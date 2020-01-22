%% Prepare logged ADC values
% 5 movement classes (Fist, Rest, Open-Hand, Wave-In, Wave-Out)
% 5 trials each
% Movements last 5s and 5s of rest in between (movement-rest-movement)
clear all;

%% Load data
file_path = 'C:\Users\JLP\Desktop\DataSync_JP'; 
str = sprintf('%s\\Myoware_500Hz_R_5g_4.csv', file_path); % File path to logged data
data = load(str);  % ADC values logged from Tiva Board

% Recording session started at rest for 10s
% Remove first 5,000 data points
% Remove extra data points at the start until it matches expected number of data points

% 500 Hz
tempor = round((size(data,1)-125000));
data = data((tempor+1:end),:); % We should have (5*5*(5s+5s)*500Hz) = 125,000 data points

%% Convert
%emg = (data*3.3)/4095;         % Map ADC to (0V) - (+3.3V)
%emg = ((data*3.3)/4095)-1.65;  % Map ADC to (-1.65) - (+1.65V)
emg = data;                     % Don't convert ADC values to EMG values

%% Class Labels
labels = [];               % Class labels
rest = zeros(2500,1);      % Rest labels (5s*200Hz)

for i=1:5  % gesture
    for j=1:5 % trial
        temp = i*ones(2500,1);
        labels = vertcat(labels, temp); % Label each trial for each movement class
        labels = vertcat(labels, rest); % Add labels for rest in between
    end
end

%% Trim data
% Extract center-portion of data for quick testing
st = 1;
en = 2500;
cut_data = [];   % Cut 3s of data from center of each gesture trial
cut_label = [];  % Labels for the cut data

for i=1:5
    for j=1:5
        tmp = emg((st:en),:);
        mid = floor(size(tmp,1)/4);
        cut_data = vertcat(cut_data, tmp((mid:mid+1499),:)); % 1500 data points (3s of data)
        cut_label = vertcat(cut_label, i*ones(1500,1));     % Class labels
        st = st+5000;
        en = en+5000;
    end
end

%% Feature Extraction
% Extracts features from data collected from 2 MyoWare sensors
% 500Hz, 2 channels

%% Initial Parameters
fs = 500;  % Sampling frequency
wl = 100;  % 200ms window, (1s/fs) = (200ms/wl)
wi = 50;   % 50% overlap, (wl-wi) = 0.5*wl

th = 0.0005;        % Threshold for ZC and TURN calculation
shfts = [-1,0,1];   % Shift neighboring analysis window for AF calculation
st = 1+1;  
en = wl+1; 

win = floor((size(cut_data,1)-wl)/wi);     % Number of windows
Nsignals = size(cut_data, 2);              % Number of electrode channels

%% Preallocate Memory for Feature Matrices
% TD
mav = zeros(win,Nsignals);     %              ***** WARNING *****
smav = zeros(win,Nsignals);    % Scaled intensity features aren't effective
rms = zeros(win,Nsignals);     % Might be because data is scaled over only 2 sensors
srms = zeros(win,Nsignals);
wav = zeros(win,Nsignals);
z = zeros(win,Nsignals);       % ZC feature is also not effective unless you filter data
t = zeros(win,Nsignals);       % Power-line interference causes DC bias

% AF
% Only meaningful if MyoWare sensors were next to each other during data collection
af_1 = zeros(win,1);
af0 = zeros(win,1);
af1 = zeros(win,1);

% Labels
lblF = zeros(win,1);

%% Filtering
% If the Tiva board was connected to a PC, data will be influenced by power-line
% interference. Other factors of noise include movement artifacts, ambient noise, etc...
% http://www.ijsrp.org/research-paper-0517/ijsrp-p6504.pdf

%cut_data = bandpass(cut_data, [20 380], fs);   % Bandpass 20-380 Hz
%cut_data = highpass(cut_data, 60, fs);       % Highpass at 60Hz
%[b, a] = butter(3, [55/(fs/2) 65/(fs/2)], 'stop'); % 6th Order Butterworth Filter (55-65 Hz Bandstop) 
%cut_data = filter(b, a, cut_data);
%cut_data = hampel(cut_data);                   % Hampel Filter

% % Notch Filter
% f0 = 60;                % Notch frequency
% fn = fs/2;              % Nyquist frequency
% freqRatio = f0/fn;      % Ratio of notch freq. to Nyquist freq.
% 
% notchWidth = 0.1;       % Width of the notch
% 
% % Compute zeros
% notchZeros = [exp( sqrt(-1)*pi*freqRatio ), exp( -sqrt(-1)*pi*freqRatio )];
% 
% % Compute poles
% notchPoles = (1-notchWidth) * notchZeros;
% 
% b = poly(notchZeros); % Get moving average filter coefficients
% a = poly(notchPoles); % Get autoregressive filter coefficients
% 
% % Filter signal
% cut_data = filter(b ,a , cut_data);

%% Feature Extraction
for i = 1:win
    for j = 1:Nsignals
        % MAV
        %mav(i,j) = mean(abs(cut_data((st:en),j)));
        mav(i,j) = mav_calc(cut_data((st:en),j));
        
        % RMS
        rms(i,j) = mean(cut_data(st:en,j).^2).^0.5;
        
        % WAV
        dx = cut_data(st:en-1,j) - cut_data(st+1:en,j);
        %wav(i,j) = mean(abs(dx));
        wav(i,j) = wav_calc(cut_data(st:en,j));
        
        % ZC and Turns
        for k=st:en
            dj = k-st+1;
            if k~=st
                if (abs(cut_data(k,j))>th) && (abs(cut_data(k-1,j))>th)
                    if (cut_data(k,j)*cut_data(k-1,j) < 0)
                        z(i,j) = z(i,j)+1;
                    end
                end
            end
            if (k<(en-1))
                if (abs(dx(dj)) && abs(dx(dj+1)) > th)
                    if (dx(dj)*dx(dj+1) < 0)
                        t(i,j) = t(i,j)+1;
                    end
                end
            end
        end
         
    end
    
    % SMAV
    mmav = mean(mav(i,:),2);
    smav(i,:) = mav(i,:)./mmav;
    
    % SRMS
    mrms = mean(rms(i,:),2);
    srms(i,:) = rms(i,:)./mrms;
    
    % Adjacent Features
    cx = cut_data((st+1:en+1),1);
    for dt = shfts
        cx = horzcat(cx, cut_data((st+dt:en+dt),2));
    end
    cX = normalizeM(cx');
    
    af_1(i,1) = sum(abs(cX(1,:)-cX(2,:)))/wl;
    af0(i,1) = sum(abs(cX(1,:)-cX(3,:)))/wl;
    af1(i,1) = sum(abs(cX(1,:)-cX(4,:)))/wl;
    
    % Majority of data points determines the class label of each analysis window 
    lblF(i,1) = round(mean(cut_label(st:en),1));
    
    st = st+wi;
    en = en+wi;
end  

% Average ZC and Turns
z = z/wl;
t = t/(wl-1);

%% Classification with the Classification Learner App
% Manually include calculated features in the training matrix:
% (mav, smav, rms, srms, wav, z, t, ar, af_1, af0, af1)

mav = double(mav);
wav = double(wav);
CL_feats = table(mav, wav, lblF);

% Open the CL App in the APPS tab in Matlab and upload this table

%% Classification with LDA Code
% Takes the features you've selected and runs them through LDA training and
% testing with code made by Dr. Zhang. 
Full_feats = table2array(CL_feats);
num_feats = zeros(5,1);

%% Cross Validation
% Divide the data into 5 parts
for n = 1:5 % Num. of gestures
    num_feats(n) = floor(length(find(Full_feats(:,end) == n))/5);
end

for p = 1:5 % testing trials
    testing(p).trial = [];
    testing(p).labels = [];
    training(p).trial = [];
    training(p).labels = [];
end

% Store 1/5 of feats as the testing matrix and 4/5 as the testing matrix
for o = 1:5 % Num. of gestures
    temp = find(Full_feats(:,end) == o);
    s = 1;
    e = num_feats(o);
    for p = 1:5 % Num. of trials
        temp_feat = Full_feats(temp(1):temp(end), 1:end-1);
        
        testing(p).trial = vertcat(testing(p).trial, temp_feat(s:e, 1:end)); % Take 1/5 as testing data
        testing(p).labels = vertcat(testing(p).labels, o*ones((e-s+1), 1));
        
        temp_feat(s:e,:) = [];   % Remove test matrix from temp full matrix
        training(p).trial = vertcat(training(p).trial, temp_feat);  % Store the remainder as a training matrix
        training(p).labels = vertcat(training(p).labels, o*ones(size(temp_feat,1), 1));
        
        s = s + num_feats(o);
        e = e + num_feats(o);
        
        if (e>size(temp,1))
            e = size(temp, 1);
        end
        
    end
end

% Transpose data in order to match syntax in the LDA train and test functions 
for p = 1:5
    testing(p).trial = testing(p).trial';
    testing(p).labels = testing(p).labels';
    training(p).trial = training(p).trial';
    training(p).labels = training(p).labels';
end

%% Results
CA = zeros(5, 1);

for q = 1:5
    % LDA Training
    train1 = training(q).trial;
    labels1 = training(q).labels;
    win1 = size(training(q).labels, 2);
    nGes = 5;
    [Wg, Cg] = LDA_train(train1, labels1, win1, nGes);
    
    % LDA Testing
    num_test = size(testing(q).trial, 2);
    prediction = zeros(1, num_test);
    
    for n=1:num_test
        %prediction(1,n) = LDA_test(testing(q).trial(:,n), Wg, Cg);
        prediction(1,n) = LDA_Predict(testing(q).trial(:,n), Wg, Cg);
    end
    
    % Count the number of correct predictions
    num = size(prediction,2);
    counter = 0;
    
    for k=1:num
        if(prediction(1,k) == testing(q).labels(1,k))
            counter = counter + 1;
        end
    end
    
    % Store classification accuracy from each cross validation fold
    CA(q) = (counter/num)*100;
end

% Average the results for 5f1 cross validation
ACA = mean(CA)

%% Save Wg and Cg in a Text File
% Training Using Full Features
[Wg, Cg] = LDA_train(Full_feats(:, 1:end-1)', lblF', size(Full_feats,1), size(unique(lblF),1));

fid = fopen('Weights\\Wg.txt','wt');

fprintf(fid, '{');
for i = 1:numel(Wg)
    fprintf(fid,'%4.5f',Wg(i));
    if (i ~= numel(Wg))
        fprintf(fid,', ');
    end
end
fprintf(fid, '}');
fclose(fid);

fim = fopen('Weights\\Cg.txt','wt');
fprintf(fim, '{');
for i = 1:size(Cg,2)
    fprintf(fim,'%4.5f',Cg(i));
    if (i ~= size(Cg,2))
        fprintf(fim, ', ');
    end
end
fprintf(fim, '}');
fclose(fim);
