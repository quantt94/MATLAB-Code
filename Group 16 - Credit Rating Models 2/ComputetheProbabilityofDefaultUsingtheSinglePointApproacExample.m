%% Compute the Probability of Default Using the Single-Point Approach to the Merton Model  

%% 
% Load the data from |MertonData.mat|. 
% load MertonData.mat
% Equity    = MertonData.Equity;
% EquityVol = MertonData.EquityVol;
% Liability = MertonData.Liability;
% Drift     = MertonData.Drift;
% Rate      = MertonData.Rate;
% MertonData  

%% 
% Compute the default probability using the single-point approach to the
% Merton model. 
clc;
clear;
Equity = 194500;
EquityVol = 0.301652;
Liability = 154191;
Rate = 0.0233;

[PD,DD,A,Sa] = mertonmodel(Equity,EquityVol,Liability,Rate);
fprintf('PD of Starbucks is: %s\n', PD);

