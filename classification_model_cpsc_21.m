clc; clear all;
load('RBF_SVM_cpsc_21.mat');
load('train_feature.mat');
load('train_label.mat');
load('test_feature.mat');
load('test_label.mat');
load("pr_test.mat");
load('indicator_test.mat');

model = fitcsvm(xtrain, ytrain, 'Standardize',true,'KernelFunction','RBF',...
    'KernelScale',1.5);
ypredict = predict(model, xtest);

% P/R and dRR based ruling
for i = 1:length(ytest)
  if ypredict(i)==0
      ypredict(i) =0;
  else
      if indicator(i)==1
          ypredict(i) = 0;
      else
          ypredict(i) = 1;
      end
  end
end

cm =confusionmat(ytest, ypredict);

% Extract values from the confusion matrix
TP = cm(1, 1);
FN = cm(1, 2);
TN = cm(2, 2);
FP = cm(2, 1);

% Accuracy
accuracy = (TP + TN) / sum(cm(:))*100;

% Sensitivity (True Positive Rate or Recall)
sensitivity = TP / (TP + FN)*100;

% Specificity (True Negative Rate)
specificity = TN / (TN + FP)*100;

% Precision (Positive Predictive Value)
precision = TP / (TP + FP)*100;

%F1 Score
f1_score = (2*precision*sensitivity)/(precision+sensitivity);
% Matthews Correlation Coefficient (MCC)
mcc = (TP * TN - FP * FN) / sqrt((TP + FP) * (TP + FN) * (TN + FP) * (TN + FN));


result =[accuracy sensitivity specificity precision f1_score mcc]