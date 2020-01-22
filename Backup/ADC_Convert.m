%% Converts logged ADC values into equivalent voltage values
% 5 movement classes (Fist, Rest, Open-Hand, Wave-In, Wave-Out)
% 5 trials each
% Movements last 5s and 5s of rest in between (movement-rest-movement)
clear all;

%% Load data
file_path = 'C:\Users\JLP\Desktop\DataSync_JP'; 
str = sprintf('%s\\Myoware_500Hz_R_5g_4.csv', file_path); % File path to logged data
data = load(str);  % ADC values logged from Tiva Board

% Started at rest for 10s, Remove first 5,000 data points
% First data point is also useless
% Remove extra data points at the start until it matches expected number of data points

% 500 Hz
tempor = round((size(data,1)-125000));
data = data((tempor+1:end),:); % We should have (5*5*(5s+5s)*500Hz) = 125,000
%% Convert
%emg = (data*3.3)/4095;        % Map ADC to (0V) - (+3.3V)
%emg = ((data*3.3)/4095)-1.65;  % Map ADC to (-1.65) - (+1.65V)
emg = data;                    % Don't convert ADC values to EMG values
%% Class Labels
labels = [];               % Class labels
rest = zeros(2500,1);      % Rest labels (5s*1000Hz) or (5s*200Hz)

for i=1:5  % gesture
    for j=1:5 % trial
        temp = i*ones(2500,1);
        labels = vertcat(labels, temp); % Label each trial for each movement class
        labels = vertcat(labels, rest); % Add labels for rest in between
    end
end

%% Extract center portion
% Used for quick testing in the next Matlab script
st = 1;
en = 2500;
cut_data = [];   % Cut 1s of data from center of each gesture trial
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

save('Data\\myoware_500Hz_R_5g_4.mat','data','emg','labels','cut_data','cut_label');
