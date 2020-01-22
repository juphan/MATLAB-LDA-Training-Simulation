%%-------LDA training procedure---------%%
% Inputs: 
%     TrainFeature: the feature matrix of the training data (feature vector dimension x # of windows)
%     TrainClass: the class label for the feature matrix (feature vector dimension x # of windows)
%     wins: total number of windows in the training dataset (number of rows in TrainFeature)
%     classes: number of classes (feature vector dimension = classes * 4)
% Outputs:
%     Wg: parameter of the LDA classifier (feature vector dimension x # of classes)
%     Cg: parameter of the LDA classifier (1 x # of classes)
%     
%%% By Xiaorong Zhang, 7/28/2014  %%%%%%%%

function [Wg, Cg]=LDA_train(TrainFeature,TrainClass,wins,classes)


%%-- Compute the means and the common covariance matrix --%%
% C: covariance matrix (feature vector dimension x feature vector dimension)
% Mi: mean vectors (feature vector dimension x # of classes)
feaure_dim=size(TrainFeature,1); % feature vector dimension
C = zeros(feaure_dim,feaure_dim); % initialize the covariance matrix
for l = 1:classes
	idx = find(TrainClass==l);
	Mi(:,l) = mean(TrainFeature(:,idx)')'; % for the data in each class, calculate a mean vector
	C = C + cov((TrainFeature(:,idx)-Mi(:,l)*ones(1,length(idx)))'); % calcuate the common covariance matrix
end
C = C./classes;
F=chol(C,'lower'); %cholesky decomposition

%%-- Compute the LDA weights --%%
for m = 1:classes
    FwL(:,m)=F\Mi(:,m);
    Wg(:,m)=F'\FwL(:,m);
    Cg(:,m) = -1/2*Mi(:,m)'*Wg(:,m);
end


 