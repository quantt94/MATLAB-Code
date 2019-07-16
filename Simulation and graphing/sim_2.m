clc;
clear all;
format long

% Initial capital
v0 = 1000;

% Number of scenarios
Ns = 5000;

% Generate Normal random variables
% normrnd(mean, stdev, numRows, numColumns)
r_speriod30 = normrnd(0.0879, 0.1465, Ns, 30);

% Distribution of value at the end of year 30
v30 = v0 * prod(1 + r_speriod30, 2);

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
[frequencyCounts, binLocations] = hist(v30, 100);
bar(binLocations, frequencyCounts);
xlabel('Value after 30 years')
ylabel('Frequency')

% Plot simulated paths over time
figure(2)
set(gcf, 'color', 'white');
time = 0:1:30;
v_t = v0*ones(Ns,1);
for(t=1:30)
    v_t = [v_t v0 * prod(1 + r_speriod30(:,1:t), 2)];
end
plot(time,v_t,'Linewidth',2);
xlabel('Time');
ylabel('Value');
title('Simulated Value Paths', 'FontWeight', 'bold');
grid on