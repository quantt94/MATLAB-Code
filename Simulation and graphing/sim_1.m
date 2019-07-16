clc;
clear all;
format long

% Initial capital
v0 = 1000;

% Number of scenarios
Ns = 100;

% Generate Normal random variables
% normrnd(mean, stdev, numRows, numColumns)
r01 = normrnd(0.0879, 0.1465, Ns, 1);

% Distribution of value at the end of year 1
v1 = (1 + r01) * v0;

% Compute statistical measures from the distribution
mean(v1) % mean
std(v1)  % standard deviation
min(v1)  % min
max(v1)  % max

% Compute percentiles/quantiles
quantile(v1, 0.05) % 5th percentile
quantile(v1,[0.05 0.50 0.95]) % 5th percentile, median and 95th percentile
sortedScen = sort(v1);  % sort scenarios
mean(sortedScen(5:6))   %  5th percentile
mean(sortedScen(95:96)) % 95th percentile
mean(sortedScen(50:51)) % median
% Alternative way to compute percentiles/quantiles
sortedScen(Ns-(1-0.05)*Ns) %  5th percentile
sortedScen(Ns-(1-0.95)*Ns) % 95th percentile

% Plot a histogram of the distribution of outcomes for v1
figure(1)
set(gcf, 'color', 'white');
% hist(v1, 10);
[frequencyCounts, binLocations] = hist(v1, 10);
bar(binLocations, frequencyCounts);
xlabel('Value at time 1')
ylabel('Frequency')

% Plot simulated paths over time
figure(2)
set(gcf, 'color', 'white');
time = 0:1:1;
plot(time,[v0*ones(Ns,1) v1],'Linewidth',2);
xlabel('Time');
ylabel('Value');
title('Simulated Value Paths', 'FontWeight', 'bold');
grid on