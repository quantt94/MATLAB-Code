clc;
clear all;
format long

% Initial capital
v0 = 1000;
% Number of scenarios
Ns = 5000;

% Expected return and covariance matrix
mu = [0.0879; 0.04];
sigma = [0.1465^2, -0.1465*0.07*0.2; -0.1465*0.07*0.2, 0.07^2];

% Generate scenarios as correlated Normal random variables
stockRet = ones(Ns,1);
bondsRet = ones(Ns,1);
for iYear = 1:30
    scenarios = mvnrnd(mu, sigma, Ns);
    stockRet = stockRet .* (1 + scenarios(:,1));
    bondsRet = bondsRet .* (1 + scenarios(:,2));
end

% Compute portfolios by iterating through different combinations of weights
for w=0.2:0.2:1
    v30(:,round(w*5)) = w*v0*stockRet + (1-w)*v0*bondsRet; % distribution of value at the end of year 30
end

% Compute statistical measures from the distribution
mean(v30) % mean
std(v30)  % standard deviation
min(v30)  % min
max(v30)  % max

% Compute percentiles/quantiles
quantile(v30, 0.05) % 5th percentile
quantile(v30,[0.05 0.50 0.95]) % 5th percentile, median and 95th percentile

% Plot a histogram of the distribution of outcomes (differences)
% for v30 (Stratery 4 - Strategy 2)
figure(1)
set(gcf, 'color', 'white');
[frequencyCounts, binLocations] = hist(v30(:,4)-v30(:,2), 50);
bar(binLocations, frequencyCounts);
xlabel('Difference in value after 30 years')
ylabel('Frequency')