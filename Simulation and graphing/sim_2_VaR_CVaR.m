clc;
clear all;
format long

% Initial capital
v0 = 1000;

% Number of scenarios
N = 5000;

% Generate Normal random variables
% normrnd(mean, stdev, numRows, numColumns)
r_speriod3 = normrnd(0.0879, 0.1465, N, 3);

% Distribution of value at the end of year 3
v3 = v0 * prod(1 + r_speriod3, 2);
% Distribution of return over 3 years
r3 = prod(1 + r_speriod3, 2) - 1;

% Losses (value and return)
loss_v3 = sort(-(v3-v0));
loss_r3 = sort(-r3);

% Quantile levels (90%, 95%, 99%, 99.9%)
alphas = [0.90 0.95 0.99 0.999];

% Compute VaR and CVaR
fprintf('Loss in value after 3 years:\n')
for(q=1:length(alphas))
    alf = alphas(q);
    VaRv(q) = loss_v3(ceil(N*alf));
    VaRr(q) = loss_r3(ceil(N*alf));
    CVaRv(q)  = (1/(N*(1-alf))) * ( (ceil(N*alf)-N*alf) * VaRv(q) + sum(loss_v3(ceil(N*alf)+1:N)) );
    CVaRr(q)  = (1/(N*(1-alf))) * ( (ceil(N*alf)-N*alf) * VaRr(q) + sum(loss_r3(ceil(N*alf)+1:N)) );
    fprintf('VaR %4.1f%% = $%6.2f, CVaR %4.1f%% = $%6.2f\n', 100*alf, VaRv(q), 100*alf, CVaRv(q))
end
fprintf('\nLoss return over 3 years:\n')
for(q=1:length(alphas))
    fprintf('VaR %4.1f%% = %6.2f%%, CVaR %4.1f%% = %6.2f%%\n', 100*alphas(q), 100*VaRr(q), 100*alphas(q), 100*CVaRr(q))
end

% Plot a histogram of the distribution of losses in value after 3 years
figure(1)
set(gcf, 'color', 'white');
[frequencyCounts, binLocations] = hist(loss_v3, 100);
bar(binLocations, frequencyCounts);
hold on;
for(q=1:length(alphas))
    line([VaRv(q) VaRv(q)], [0 max(frequencyCounts)], 'Color', 'r', 'LineWidth', 1, 'LineStyle', '--')
    hold on;
end
hold off;
xlabel('Loss in portfolio value after 3 years')
ylabel('Frequency')