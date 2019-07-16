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

% Generate correlated Normal random variables
stockRet = ones(Ns,1);
bondsRet = ones(Ns,1);
v_t = v0 * ones(Ns,1);
for iYear = 1:30
    scenarios = mvnrnd(mu, sigma, Ns);
    stockRet = stockRet .* (1 + scenarios(:,1));
    bondsRet = bondsRet .* (1 + scenarios(:,2));
    v_t = [v_t 0.5*v0*stockRet + 0.5*v0*bondsRet];
end

% Distribution of value at the end of year 30
v30 = 0.5*v0*stockRet + 0.5*v0*bondsRet;

% Compute statistical measures from the distribution
mean(v30) % mean
std(v30)  % standard deviation
min(v30)  % min
max(v30)  % max

% Compute percentiles/quantiles
quantile(v30, 0.05) % 5th percentile
quantile(v30,[0.05 0.50 0.95]) % 5th percentile, median and 95th percentile
sortedScen = sort(v30);  % sort scenarios
mean(sortedScen(5:6))   %  5th percentile
mean(sortedScen(95:96)) % 95th percentile
mean(sortedScen(50:51)) % median
% Alternative way to compute percentiles/quantiles
sortedScen(Ns-(1-0.05)*Ns) %  5th percentile
sortedScen(Ns-(1-0.95)*Ns) % 95th percentile

% Plot a histogram of the distribution of outcomes for v30
figure(1)
set(gcf, 'color', 'white');
% hist(v30, 10);
[frequencyCounts, binLocations] = hist(v30, 50);
bar(binLocations, frequencyCounts);
xlabel('Value after 30 years')
ylabel('Frequency')

% Plot simulated paths over time
figure(2)
set(gcf, 'color', 'white');
time = 0:1:30;
plot(time,v_t,'Linewidth',2);
xlabel('Time');
ylabel('Value');
title('Simulated Value Paths', 'FontWeight', 'bold');
grid on