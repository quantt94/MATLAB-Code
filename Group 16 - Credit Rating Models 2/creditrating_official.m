clc;
clear;
%% Import the data
[~, ~, raw] = xlsread('/Users/JunJun/Documents/MATLAB/training.xls','Sheet1');
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};
stringVectors = string(raw(:,[1,2,3,16]));
stringVectors(ismissing(stringVectors)) = '';
raw = raw(:,[4,5,6,7,8,9,10,11,12,13,14,15]);

%% Exclude rows with non-numeric cells
I = ~all(cellfun(@(x) (isnumeric(x) || islogical(x)) && ~isnan(x),raw),2); % Find rows with non-numeric cells
raw(I,:) = [];
stringVectors(I,:) = [];

%% Create output variable
data = reshape([raw{:}],size(raw));

%% Create table
dataset1 = table;

%% Allocate imported array to column variable names
dataset1.TICKER = stringVectors(:,1);
dataset1.SECURITY_NAME = stringVectors(:,2);
dataset1.RTG_SP_LT_LC_ISSUER_CREDIT = categorical(stringVectors(:,3));
dataset1.BS_TOT_ASSET = data(:,1);
dataset1.BS_TOT_LIAB2 = data(:,2);
dataset1.TOTAL_EQUITY = data(:,3);
dataset1.HISTORICAL_MARKET_CAP = data(:,4);
dataset1.WORKING_CAPITAL = data(:,5);
dataset1.BS_PURE_RETAINED_EARNINGS = data(:,6);
dataset1.RETAIN_EARN_TO_TOT_ASSET = data(:,7);
dataset1.EBIT = data(:,8);
dataset1.SALES_REV_TURN = data(:,9);
dataset1.ALTMAN_Z_SCORE = data(:,10);
dataset1.BS_CUR_ASSET_REPORT = data(:,11);
dataset1.BS_CUR_LIAB = data(:,12);
dataset1.INDUSTRY_SECTOR = categorical(stringVectors(:,4));

%% Clear temporary variables
clearvars data raw stringVectors I;
%% CALCULATE FACTORS 
WC_TA = dataset1.WORKING_CAPITAL./dataset1.BS_TOT_ASSET;

RE_TA = dataset1.RETAIN_EARN_TO_TOT_ASSET;

EBIT_TA = dataset1.EBIT./dataset1.BS_TOT_ASSET;

MVE_BVTD = dataset1.HISTORICAL_MARKET_CAP./dataset1.BS_TOT_LIAB2;

S_TA = dataset1.SALES_REV_TURN./dataset1.BS_TOT_ASSET;

Industry = grp2idx(dataset1.INDUSTRY_SECTOR);

X = [WC_TA, ... 
     RE_TA, ...
     EBIT_TA, ... 
     MVE_BVTD, ...
     S_TA, ...
     Industry];
 
% Set aside one observation to test the model later 
    rng('default');
    ind = randi([1,length(X)-3],1,length(X)-3);
    X_train = X(ind,:);

Y = ordinal(dataset1.RTG_SP_LT_LC_ISSUER_CREDIT(ind),[],{'AAA' 'AA+' 'AA' 'AA-'...
    'A+' 'A' 'A-' 'BBB+' 'BBB' 'BBB-' 'BB+' 'BB' 'BB-' 'B+' 'B' 'B-' 'CCC+'...
    'CCC' 'CCC-' 'CC' 'C' 'RD' 'SD' 'D'});

%%
% We will also create a numerical equivalent of |Y| (where AAA = 1, AA = 2,
% etc.) for the benefit of some of the simpler regression functions.
Y_num = double(Y);

%% Choose number of trees
figure;
leaf = [1];
nTrees = 70;
rng(1128,'twister');
savedRng = rng; % save the current RNG settings

color = 'bgr';
for ii = 1:length(leaf)
   % Reinitialize the random number generator, so that the
   % random samples are the same for each leaf size
   rng(savedRng);
   % Create a bagged decision tree for each leaf size and plot out-of-bag
   % error 'oobError'
   b = TreeBagger(nTrees,X_train,Y,'OOBPred','on',...
                             'CategoricalPredictors',6,...
                             'MinLeaf',leaf(ii));
   plot(b.oobError,color(ii));
   hold on;
end
xlabel('Number of grown trees');
ylabel('Out-of-bag classification error');
title('Classification Error for Different Leaf Sizes');
hold off;

%%
% The errors are comparable for the three leaf-size options. We will
% therefore work with a leaf size of 10, because it results in leaner trees
% and more efficient computations.
%
% Note that we did not have to split the data into _training_ and _test_
% subsets. This is done internally, it is implicit in the sampling
% procedure that underlies the method. At each bootstrap iteration, the
% bootstrap replica is the training set, and any customers left out
% ("out-of-bag") are used as test points to estimate the out-of-bag
% classification error reported above.
%
% Next, we want to find out whether all the features are important for the
% accuracy of our classifier. We do this by turning on the _feature
% importance_ measure (|oobvarimp|), and plot the results to visually find
% the most important features. We also try a larger number of trees now,
% and store the classification error, for further comparisons below.
figure;
nTrees = 35;
leaf = 1;
rng(1128,'twister');
b = TreeBagger(nTrees,X_train,Y,'oobvarimp','on','cat',6,'minleaf',leaf);

figure(2);
bar(b.OOBPermutedVarDeltaError);
xlabel('Feature number');
ylabel('Out-of-bag feature importance');
title('Feature importance results');

oobErrorFullX = b.oobError;
%%
Y_tb = b.predict(X_train);
Y_tb = ordinal(Y_tb,[],{'AAA' 'AA+' 'AA' 'AA-'...
    'A+' 'A' 'A-' 'BBB+' 'BBB' 'BBB-' 'BB+' 'BB' 'BB-' 'B+' 'B' 'B-' 'CCC+'...
    'CCC' 'CCC-' 'CC' 'C' 'RD' 'SD' 'D'});
%%
starbuck = X(364,:);
actual_credit = dataset1.RTG_SP_LT_LC_ISSUER_CREDIT(364);

%%
% To predit the credit rating for this new data, we call the |predict|
% method on the classifier. The method returns two arguments, the predicted
% class and the classification score. We certainly want to get both output
% arguments, since the classification scores contain information on how
% certain the predicted ratings seem to be.

[predClass,classifScore] = b.predict(starbuck);
%%
% At this point, we can create a report. Here we only display a small
% report for the first three customers on the screen, for illustration
% purposes, but a more detailed report could be written to a file as well.

   disp(dataset1.SECURITY_NAME(364));
   fprintf('     Industry: %s\n',char(dataset1.INDUSTRY_SECTOR(364)));
   fprintf('     Predicted Rating: %s\n',char(predClass));
   fprintf('     Actual Rating: %s\n', char(actual_credit));
