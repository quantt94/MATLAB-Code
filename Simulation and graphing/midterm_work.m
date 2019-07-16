clc;
clear all;
format long

% CSV file with price data
input_file_prices  = 'Daily_closing_prices.csv';

% Read daily prices
if(exist(input_file_prices,'file'))
  fprintf('\nReading daily prices datafile - %s\n', input_file_prices)
  fid = fopen(input_file_prices);
     % Read instrument tickers
     hheader  = textscan(fid, '%s', 1, 'delimiter', '\n');
     headers = textscan(char(hheader{:}), '%q', 'delimiter', ',');
     tickers = headers{1}(2:end);
     % Read time periods
     vheader = textscan(fid, '%[^,]%*[^\n]');
     dates = vheader{1}(1:end);
  fclose(fid);
  data_prices = dlmread(input_file_prices, ',', 1, 1);
else
  error('Daily prices datafile does not exist')
end

% Convert dates into array [year month day]
format_date = 'mm/dd/yyyy';
dates_array = datevec(dates, format_date);
dates_array = dates_array(:,1:3);

% Remove datapoints for year 2014
day_ind_start0 = 1;
day_ind_end0 = length(find(dates_array(:,1)==2014));
data_prices = data_prices(day_ind_end0+1:end,:);
dates_array = dates_array(day_ind_end0+1:end,:);
dates = dates(day_ind_end0+1:end,:);

% Compute means and covariances for Question 2
day_ind_start = 1;
day_ind_end = 39;
his_returns = data_prices(day_ind_start+1:day_ind_end,:) ./ data_prices(day_ind_start:day_ind_end-1,:) - 1;
mu = mean(his_returns)';  % Expected returns for Question 2
Q = cov(his_returns);     % Covariances for Question 2


%% Question 1 - Part I

% Specify quantile level for VaR/CVaR
alf = 0.95;

% Positions in the portfolio
positions = [100 0 0 0 0 0 0 0 200 500 0 0 0 0 0 0 0 0 0 0]';

% Number of assets in universe
Na = size(data_prices,2);

% Number of historical scenarios
Ns = size(data_prices,1);

%%%%% Insert your code here 

% Calculate Daily Portfolio value
port_value = data_prices*positions;

% Number of Historical scenarios for 1-day horizon
N1d = length(port_value) - 1;

% Number of historical scenarios for 10-day horizon

N10d = length(port_value) - 9;

% Calculate VaR and CVaR for 1 day horizon 

    % create the vector to store P&L
    pl_1d = [];
    
    % create the loop to caculate P&L for 1 day horizon
    for i = 1:N1d;
        pl_1d(i) = port_value(i+1) - port_value(i);
    end
    
    % sort the loss
    loss_1d = sort(-pl_1d);
    
    % compute historical 1-day VaR and CVaR from the data
    VaR1  = loss_1d(ceil(N1d*alf));
    CVaR1 = (1/(N1d*(1-alf))) * ( (ceil(N1d*alf)-N1d*alf) * VaR1 + sum(loss_1d(ceil(N1d*alf)+1:N1d)) );
    
    
    cur_value = data_prices(504,:).*positions';
    index = find(cur_value > 0);
    new_data_prices = data_prices(:,index);
    new_positions = positions(index);
        % calculate loss of each stocks and take the mean and variance
        % covariance 
        loss_1d_normal = [];
        mean_loss_1d_normal = [];
        for i = 1:3
            for j = 1:N1d
                loss_1d_normal(j,i) = -(new_data_prices(j+1,i) - new_data_prices(j,i));
            end
        end
        mean_loss_1d_normal = mean(loss_1d_normal);
        var_loss_1d_normal = cov(loss_1d_normal);
        
        % calculate mean of portfolio loss and variance of portfolio loss
        mean_loss_1d = mean_loss_1d_normal * new_positions;
        var_loss_1d = new_positions' * var_loss_1d_normal * new_positions;
        stdev_loss_1d = sqrt(var_loss_1d);
    
        
        
    % compute normal 1-day VaR and CVaR from the data
    VaR1n = mean_loss_1d + norminv(alf,0,1)*stdev_loss_1d;
    CVaR1n = mean_loss_1d + (normpdf(norminv(alf,0,1))/(1-alf))*stdev_loss_1d;
    
% Calculate VaR and CVaR for 10 day horizon

    % create the vector to store P&L
    pl_10d = zeros(N10d,1);
    
    % create the loop to caculate P&L for 10 day horizon
    for i = 1:N10d;
        pl_10d(i) = port_value(i+9) - port_value(i);
    end
    
    % sort the loss
    loss_10d = sort(-pl_10d);
    
    % compute historical 1-day VaR and CVaR from the data
    VaR10  = loss_10d(ceil(N10d*alf));
    CVaR10 = (1/(N10d*(1-alf))) * ( (ceil(N10d*alf)-N10d*alf) * VaR10 + sum(loss_10d(ceil(N10d*alf)+1:N10d)) );
    
    % calculate loss of each stocks and take the mean and variance
        % covariance 
        loss_10d_normal = [];
        mean_loss_10d_normal = [];
        for i = 1:3
            for j = 1:N10d
                loss_10d_normal(j,i) = -(new_data_prices(j+9,i) - new_data_prices(j,i));
            end
        end
        mean_loss_10d_normal = mean(loss_10d_normal);
        var_loss_10d_normal = cov(loss_10d_normal);
        
        % calculate mean of portfolio loss and variance of portfolio loss
        mean_loss_10d = mean_loss_10d_normal * new_positions;
        var_loss_10d = new_positions' * var_loss_10d_normal * new_positions;
        stdev_loss_10d = sqrt(var_loss_10d);
        
    % compute normal 10-day VaR and CVaR from the data
    VaR10n = mean_loss_10d + norminv(alf,0,1)*stdev_loss_10d;
    CVaR10n = mean_loss_10d + (normpdf(norminv(alf,0,1))/(1-alf))*stdev_loss_10d;
fprintf('QUESTION 1 - PART 1 \n')   
fprintf('    Historical 1-day VaR %4.1f%% = $%6.2f,   Historical 1-day CVaR %4.1f%% = $%6.2f\n', 100*alf, VaR1, 100*alf, CVaR1)
fprintf('        Normal 1-day VaR %4.1f%% = $%6.2f,       Normal 1-day CVaR %4.1f%% = $%6.2f\n', 100*alf, VaR1n, 100*alf, CVaR1n)
fprintf ('\n')
fprintf('    Historical 10-day VaR %4.1f%% = $%6.2f,   Historical 10-day CVaR %4.1f%% = $%6.2f\n', 100*alf, VaR10, 100*alf, CVaR10)
fprintf('        Normal 10-day VaR %4.1f%% = $%6.2f,       Normal 10-day CVaR %4.1f%% = $%6.2f\n', 100*alf, VaR10n, 100*alf, CVaR10n)
fprintf('\n')

% Plot a histogram of the distribution of losses in portfolio value for 1 day 
% figure(1)
    figure(1)

% set the figure size to be bigger
    set(gcf, 'color','white','Units', 'Normalized', 'OuterPosition', [0, 0.03, 0.7, 0.8]);
    set(gca, 'FontName', 'Calibri');

    [frequencyCounts, binLocations] = hist(loss_1d, 100);
    bar(binLocations, frequencyCounts,'b');
hold on;
    line([VaR1 VaR1], [0 max(frequencyCounts)/2], 'Color', [0.9100 0.4100 0.1700], 'LineWidth', 1.5, 'LineStyle', '-');
hold on;
    normf = ( 1/(std(loss_1d)*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(loss_1d))/std(loss_1d)).^2 );
    normf = normf * sum(frequencyCounts)/sum(normf);
    plot(binLocations, normf, 'r', 'LineWidth', 3);
hold on;
    line([VaR1n VaR1n], [0 max(frequencyCounts)/2], 'Color', [0 0.5 0], 'LineWidth', 1.5, 'LineStyle', '-');
hold on;
    line([CVaR1 CVaR1], [0 max(frequencyCounts)/2], 'Color', [0.9100 0.4100 0.1700], 'LineWidth', 1.5, 'LineStyle', '--');
hold on;
    line([CVaR1n CVaR1n], [0 max(frequencyCounts)/2], 'Color', [0 0.5 0], 'LineWidth', 1.5, 'LineStyle', '--');
hold off;
    txt = '\leftarrow VaR Historical';
    text(1.05*VaR1,max(frequencyCounts)/2,txt,'HorizontalAlignment','left')
    txt = 'VaR Normal \rightarrow';
    text(0.95*VaR1n,max(frequencyCounts)/2,txt,'HorizontalAlignment','right')
    txt = '\leftarrow CVaR Historical';
    text(1.5*VaR1,max(frequencyCounts)/2.5,txt,'HorizontalAlignment','left')
    txt = 'CVaR Normal \rightarrow';
    text(1.2*VaR1n,max(frequencyCounts)/2.5,txt,'HorizontalAlignment','right')
    
%     text(0.8*VaR1, max(frequencyCounts)/1.9, 'VaR1', 'FontSize', 12, 'FontWeight','bold')
%     text(0.98*VaR1n, max(frequencyCounts)/1.9, 'VaR1n','FontSize', 12,'FontWeight','bold')

title('95% VaR AND CVaR IN DOLLAR VALUE OF DAILY PORTFOLIO VALUE','FontSize', 18);
    xlabel('1-day Loss in Dollar Value','FontWeight','bold','FontSize',14);
    ylabel('Frequency', 'FontWeight','bold','FontSize',14);
lgd = legend('Historical Distribution','Historical VaR','Normal Distribution','Normal VaR','Historical CVaR'...
    ,'Normal CVaR','location','northwest');
    lgd.Title.String = 'METHOD';
    lgd.FontSize = 14;

% Plot a histogram of the distribution of losses in portfolio value for 10 days
% figure(2)
    figure(2)

% set the figure size to be bigger
    set(gcf, 'color','white','Units', 'Normalized', 'OuterPosition', [0, 0.03, 0.7, 0.8]);
    set(gca, 'FontName', 'Calibri');

    [frequencyCounts, binLocations] = hist(loss_10d, 100);
    bar(binLocations, frequencyCounts,'b');
hold on;
    line([VaR10 VaR10], [0 max(frequencyCounts)/2], 'Color', [0.9100 0.4100 0.1700], 'LineWidth', 1.5, 'LineStyle', '-');
hold on;
    normf = ( 1/(std(loss_10d)*sqrt(2*pi)) ) * exp( -0.5*((binLocations-mean(loss_10d))/std(loss_10d)).^2 );
    normf = normf * sum(frequencyCounts)/sum(normf);
    plot(binLocations, normf, 'r', 'LineWidth', 3);
hold on;
    line([VaR10n VaR10n], [0 max(frequencyCounts)/2], 'Color', [0 0.5 0], 'LineWidth', 1.5, 'LineStyle', '-');
hold on;
    line([CVaR10 CVaR10], [0 max(frequencyCounts)/2], 'Color', [0.9100 0.4100 0.1700], 'LineWidth', 1.5, 'LineStyle', '--');
hold on;
    line([CVaR10n CVaR10n], [0 max(frequencyCounts)/2], 'Color', [0 0.5 0], 'LineWidth', 1.5, 'LineStyle', '--');
hold off;
    txt = '\leftarrow VaR Historical';
    text(1.02*VaR10,max(frequencyCounts)/2,txt,'HorizontalAlignment','left')
    txt = 'VaR Normal \rightarrow';
    text(0.98*VaR10n,max(frequencyCounts)/2,txt,'HorizontalAlignment','right')
    txt = '\leftarrow CVaR Historical';
    text(1.25*VaR10,max(frequencyCounts)/2.5,txt,'HorizontalAlignment','left')
    txt = 'CVaR Normal \rightarrow';
    text(1.23*VaR10n,max(frequencyCounts)/2.5,txt,'HorizontalAlignment','right')
    
title('95% VaR AND CVaR IN DOLLAR VALUE OF DAILY PORTFOLIO VALUE','FontSize', 18);
    xlabel('10-day Loss in Dollar Value','FontWeight','bold','FontSize',14);
    ylabel('Frequency', 'FontWeight','bold','FontSize',14);
lgd = legend('Historical Distribution','Historical VaR','Normal Distribution'...
    ,'Normal VaR','Historical CVaR','Normal CVaR','location','northwest');
    lgd.Title.String = 'METHOD';
    lgd.FontSize = 14;

%% Question 1 - Part 2

% Calculate Daily Portfolio value assuming only MSFT 
MSFT = data_prices(:,1).*positions(1);

% Calculate VaR and CVaR for 1 day horizon 

    % create the vector to store P&L
    pl_MSFT = zeros(N1d,1);
    
    % create the loop to caculate P&L for 1 day horizon
    for i = 1:N1d;
        pl_MSFT(i) = MSFT(i) - MSFT(i+1);
    end
    
    % sort the loss
    loss_MSFT = sort(pl_MSFT);
    
    % compute historical 1-day VaR and CVaR from the data
    VaR_MSFT  = loss_MSFT(ceil(N1d*alf));
    CVaR_MSFT = (1/(N1d*(1-alf))) * ( (ceil(N1d*alf)-N1d*alf) * VaR_MSFT + sum(loss_MSFT(ceil(N1d*alf)+1:N1d)) );
    % compute normal 1-day VaR and CVaR from the data
    VaRn_MSFT = mean(loss_MSFT) + norminv(alf,0,1)*std(loss_MSFT);
    CVaRn_MSFT = mean(loss_MSFT) + (normpdf(norminv(alf,0,1))/(1-alf))*std(loss_MSFT);
    
% Calculate Daily Portfolio value assuming only APPL 
APPL = data_prices(:,9).*positions(9);

% Calculate VaR and CVaR for 1 day horizon 

    % create the vector to store P&L
    pl_APPL = zeros(N1d,1);
    
    % create the loop to caculate P&L for 1 day horizon
    for i = 1:N1d;
        pl_APPL(i) = APPL(i) - APPL(i+1);
    end
    
    % sort the loss
    loss_APPL = sort(pl_APPL);
    
    % compute historical 1-day VaR and CVaR from the data
    VaR_APPL  = loss_APPL(ceil(N1d*alf));
    CVaR_APPL = (1/(N1d*(1-alf))) * ( (ceil(N1d*alf)-N1d*alf) * VaR_APPL + sum(loss_APPL(ceil(N1d*alf)+1:N1d)) );
    % compute normal 1-day VaR and CVaR from the data
    VaRn_APPL = mean(loss_APPL) + norminv(alf,0,1)*std(loss_APPL);
    CVaRn_APPL = mean(loss_APPL) + (normpdf(norminv(alf,0,1))/(1-alf))*std(loss_APPL);
    
% Calculate Daily Portfolio value assuming only IBM 
IBM = data_prices(:,10).*positions(10);

% Calculate VaR and CVaR for 1 day horizon 

    % create the vector to store P&L
    pl_IBM = zeros(N1d,1);
    
    % create the loop to caculate P&L for 1 day horizon
    for i = 1:N1d;
        pl_IBM(i) = IBM(i) - IBM(i+1);
    end
    
    % sort the loss
    loss_IBM = sort(pl_IBM);
    
    % compute historical 1-day VaR and CVaR from the data
    VaR_IBM  = loss_IBM(ceil(N1d*alf));
    CVaR_IBM = (1/(N1d*(1-alf))) * ( (ceil(N1d*alf)-N1d*alf) * VaR_IBM + sum(loss_IBM(ceil(N1d*alf)+1:N1d)) );
    % compute normal 1-day VaR and CVaR from the data
    VaRn_IBM = mean(loss_IBM) + norminv(alf,0,1)*std(loss_IBM);
    CVaRn_IBM = mean(loss_IBM) + (normpdf(norminv(alf,0,1))/(1-alf))*std(loss_IBM);

fprintf('QUESTION 1 - PART 2 \n')
fprintf('    Historical 1-day VaR of MSFT %4.1f%% = $%6.2f,   Historical 1-day CVaR of MSFT %4.1f%% = $%6.2f\n', 100*alf, VaR_MSFT, 100*alf, CVaR_MSFT)
fprintf('        Normal 1-day VaR of MSFT %4.1f%% = $%6.2f,       Normal 1-day CVaR of MSFT %4.1f%% = $%6.2f\n', 100*alf, VaRn_MSFT, 100*alf, CVaRn_MSFT)
fprintf('    Historical 1-day VaR of APPL %4.1f%% = $%6.2f,   Historical 1-day CVaR of APPL %4.1f%% = $%6.2f\n', 100*alf, VaR_APPL, 100*alf, CVaR_APPL)
fprintf('        Normal 1-day VaR of APPL %4.1f%% = $%6.2f,       Normal 1-day CVaR of APPL %4.1f%% = $%6.2f\n', 100*alf, VaRn_APPL, 100*alf, CVaRn_APPL)
fprintf('    Historical 1-day VaR of IBM  %4.1f%% = $%6.2f,   Historical 1-day CVaR of IBM   %4.1f%% = $%6.2f\n', 100*alf, VaR_IBM, 100*alf, CVaR_IBM)
fprintf('        Normal 1-day VaR of IBM  %4.1f%% = $%6.2f,       Normal 1-day CVaR of IBM   %4.1f%% = $%6.2f\n', 100*alf, VaRn_IBM, 100*alf, CVaRn_IBM)
fprintf(' \n')


%% Question 2 - Part 1
% Question 2

% Annual risk-free rate for years 2015-2016 is 2.5%
r_rf = 0.025;

% Initial portfolio weights
init_positions = [5000 950 2000 0 0 0 0 2000 3000 1500 0 0 0 0 0 0 1001 0 0 0]';

% this is the value of 3/2/2015 - initial weight of Assignment 1
init_value = data_prices(day_ind_end+1,:) * init_positions;
w_init = (data_prices(day_ind_end+1,:) .* init_positions')' / init_value;

% Max Sharpe Ratio portfolio weights
w_Sharpe = [ 0 0 0 0 0 0 0 0.385948690661642 0.172970428625544 0 0 0 0 0 0.003409676869715 0.260942060896445 0 0.185966939781285 0 0]';

% Equal Risk Contribution portfolio weights
w_ERC = [0.049946771209069 0.049951626261681 0.049955739901370 0.049998404150207 0.050000297368719 0.050004255546315 0.050006307026730 0.050007308995726 0.050010525832832 0.050013840015521 0.050014404492514 0.050015932843104 0.050016630302524 0.050017212457105 0.050017600497611 0.050017998351827 0.050018997074443 0.050019598350121 0.050019778113513 0.049946771209069]';

%%%%% Insert your code here 

% use historical data to estimate return and VcV matrix
% which is already given in previous code

% Optimization problem data
lb = zeros(Na,1); % lower bound with a vector of n zeros
ub = inf*ones(Na,1); % upper bound with a vector of n infinities
A  = ones(1,Na);
b  = 1;

% MINIMUM VARIANCE PORTFOLIO

    % add name to function
    cplex1 = Cplex('min_Variance');
    % Add objective function and bounds on variables to CPLEX model
    cplex1.addCols(zeros(Na,1), [], lb, ub); 
    % Add constraints to CPLEX model
    cplex1.addRows(b, A, b); 
    % Add quadratic part of objective function to CPLEX model
    cplex1.Model.Q = 2*Q;
    % concurrent algorithm
    cplex1.Param.qpmethod.Cur = 6; 
    % enable crossover
    cplex1.Param.barrier.crossover.Cur = 1; 
    % disable output to screen
    cplex1.DisplayFunc = []; 
    % solve
    cplex1.solve();

    % Display minimum variance portfolio
    w_minVar = cplex1.Solution.x;
    var_minVar = w_minVar' * Q * w_minVar;
    ret_minVar = mu' * w_minVar;
    stdev_minVar = sqrt(var_minVar);
    
    fprintf ('QUESTION 2 - PART 1 \n');
    fprintf ('    Minimum variance portfolio:\n');
    fprintf ('      Return = %f\n', ret_minVar);
    fprintf ('      Standard deviation = %f\n\n', sqrt(var_minVar));
    fprintf ('\n')

% MAXIMUM RETURN PORTFOLIO 

    cplex2 = Cplex('max_Return');
    cplex2.Model.sense = 'maximize';
    cplex2.addCols(mu, [], lb, ub);
    cplex2.addRows(b, A, b);
    cplex2.Param.lpmethod.Cur = 6; % concurrent algorithm
    cplex2.Param.barrier.crossover.Cur = 1; % enable crossover
    cplex2.DisplayFunc = []; % disable output to screen
    cplex2.solve();

    % Display maximum return portfolio
    w_maxRet = cplex2.Solution.x;
    var_maxRet = w_maxRet' * Q * w_maxRet;
    ret_maxRet = mu' * w_maxRet;
    stdev_maxRet = sqrt(var_maxRet);
    fprintf ('    Maximum return portfolio:\n');
    fprintf ('      Return = %f\n', ret_maxRet);
    fprintf ('      Standard deviation = %f\n\n', sqrt(var_maxRet));
    fprintf ('\n')

% Target returns, equally divided return-axis to 20 points between the
% return of min var port. and the return of max return port.
targetRet = linspace(ret_minVar,ret_maxRet,20);

% EFFICIENT FRONTIER OF RISKY ASSETS

    cplex3 = Cplex('Efficient_Frontier');
    cplex3.addCols(zeros(Na,1), [], lb, ub);
    cplex3.addRows(targetRet(1), mu', inf);
    cplex3.addRows(b, A, b);
    cplex3.Model.Q = 2*Q;
    cplex3.Param.qpmethod.Cur = 6; % concurrent algorithm
    cplex3.Param.barrier.crossover.Cur = 1; % enable crossover
    cplex3.DisplayFunc = []; % disable output to screen

    w_front = [];
    for i=1:length(targetRet)
        cplex3.Model.lhs(1) = targetRet(i);
        cplex3.solve();
        w_front = [w_front cplex3.Solution.x];
        var_front(i) = w_front(:,i)' * Q * w_front(:,i);
        ret_front(i) = mu' * w_front(:,i);
    end

% EQUALLY WEIGHTED PORTFOLIO 

    w_equal = (1/Na)*ones(Na,1);
    ret_equal = mu' * w_equal;
    var_equal = w_equal' * Q * w_equal;
    stdev_equal = sqrt(var_equal);
    fprintf ('    Equally weighted portfolio:\n');
    fprintf ('      Return = %f\n', ret_equal);
    fprintf ('      Standard deviation = %f\n\n', stdev_equal);
  
% INITIAL PORTFOLIO

    ret_init = mu' * w_init;
    var_init = w_init' * Q * w_init;
    stdev_init = sqrt(var_init);
    fprintf ('    Initial portfolio:\n');
    fprintf ('      Return = %f\n', ret_init);
    fprintf ('      Standard deviation = %f\n\n', stdev_init);
    
% MAXIMUM SHARPE RATIO PORTFOLIO 

    ret_Sharpe = mu' * w_Sharpe;
    var_Sharpe = w_Sharpe' * Q * w_Sharpe;
    stdev_Sharpe = sqrt(var_Sharpe);
    fprintf ('    Maximum Sharpe ratio portfolio:\n');
    fprintf ('      Return = %f\n', ret_Sharpe);
    fprintf ('      Standard deviation = %f\n\n', stdev_Sharpe);
    
% RISK FREE ASSETS
    
    ret_rf = r_rf/252;
    var_rf = 0;
    stdev_rf = 0;
    fprintf ('    Risk-free Asset:\n');
    fprintf ('      Return = %f\n', ret_rf);
    fprintf ('      Standard deviation = %f\n\n', stdev_rf);
    
% EQUAL RISK CONTRIBUTIONS 

    ret_ERC = mu' * w_ERC;
    var_ERC = w_ERC' * Q * w_ERC;
    stdev_ERC = sqrt(var_ERC);
    fprintf ('    Equal risk contributions portfolio:\n');
    fprintf ('      Return = %f\n', ret_ERC);
    fprintf ('      Standard deviation = %f\n\n', stdev_ERC);
    
% LEVERAGE EQUAL RISK CONTRIBUTIONS
    % borrowing or lending risk free
    % the line connecting risk-free assets and  passing through ERC portfolio
    y = [ret_rf ret_ERC];
    x = [stdev_rf stdev_ERC];
    slope1 = (y(2) - y(1))/(x(2) - x(1));
    intercept1 = ret_ERC - slope1 * stdev_ERC;
    
    stdev_lev = [0 stdev_ERC:0.001:stdev_maxRet];
    ret_lev = intercept1 + slope1 * stdev_lev;
    

% EFFICIENT FRONTIER WITH RISK-FREE ASSETS SHORT SALES ALLOWED BUT NOT
% TANGENT LINE EXPANDED AFTER SHARPE RATIO
    
    ret_front_rf = [ret_rf ret_Sharpe];
    stdev_front_rf = [stdev_rf stdev_Sharpe];

% Plot for Question 2, Part 1
% figure(3);
% Plot efficient frontier
figure(3);
% set the figure size to be bigger
    set(gcf, 'color','white','Units', 'Normalized', 'OuterPosition', [0, 0.03, 0.7, 0.8]);
    set(gca, 'FontName', 'Calibri');
    
    plot(sqrt(var_front), ret_front, 'k-', 'LineWidth', 3, 'color', 'b')
hold on;
    plot(sqrt(var_minVar), ret_minVar, 'r*', 'MarkerSize', 10)
hold on;
    plot(sqrt(var_maxRet), ret_maxRet, 'ks', 'MarkerSize', 10, 'MarkerFaceColor','m')
hold on;
    plot(stdev_equal, ret_equal, 'bd', 'MarkerSize', 10, 'MarkerFaceColor', 'g')
hold on;
    plot(stdev_Sharpe, ret_Sharpe, 'ko', 'MarkerSize', 10, 'MarkerFaceColor', 'y')
hold on;
    plot(stdev_rf, ret_rf, '.', 'MarkerSize', 30)
hold on;
    plot(stdev_ERC, ret_ERC, 'v', 'MarkerSize', 8, 'MarkerFaceColor', 'r', 'MarkerEdgeColor','k')
hold on;
    plot(stdev_lev, ret_lev, ':', 'color', [0.9100 0.4100 0.1700],'LineWidth', 2)
hold on;
    plot(stdev_init, ret_init, '^', 'MarkerSize', 10, 'MarkerEdgeColor', 'w', 'MarkerFaceColor','k')
hold on;
    plot(stdev_front_rf, ret_front_rf, ':', 'color', [0 0.5 0], 'LineWidth', 1.5)
    
hold off;

% display grid accord to y-axis
    grid on
    ax = gca;
    ax.GridColor = [0.1 0.1 0.1];
    ax.GridLineStyle = '--';
% name the graph, x-axis, y-axis and legend
    title('EFFICIENT FRONTIER AND DIFFERENT PORTFOLIO STRATEGIES ON MARCH 2015 ','FontSize', 18);
    xlabel('Portfolio Standard Deviation','FontWeight','bold','FontSize',14);
    ylabel('Portfolio Return', 'FontWeight','bold','FontSize',14);
    
lgd = legend('Efficient Frontier', 'Minimum Variance Portfolio', 'Maximum Return Portfolio', 'Equally Weighted Portfolio'...
    ,'Maximum Sharpe Ratio','Risk Free Asset','Equal Risk Contributions'...
    ,'Leverage Equal Risk Contributions','Initial Portfolio','Efficient Frontier with Risk Free',...
    'Location', 'NorthWest');
    lgd.Title.String = 'STRATEGY NAME';
    lgd.FontSize = 12;
    
%% Question 2 - Part 2
% RANDOMLY GENERATED 1000 PORTFOLIOS UNDER NO-SHORT-SALES RESTRICTIONS

    % create storage for mu and variance from random generator
    mu_random = zeros(1000,1);
    var_random = zeros(1000,1);
    

    % normalize the weights to sum up to one for each of the scenario
    minval = 0;
    maxval = 1;
        for i = 1:1000
            w_random = minval + (maxval - minval).*rand(20,1);
            w_random = (w_random - min(w_random))/(max(w_random) - min(w_random));
            w_random = w_random/sum(w_random);
            ret_random(i) = mu'*w_random;
            stdev_random(i) = sqrt(w_random' * Q * w_random);
        end
        
% figure(4);
figure(4);
% set the figure size to be bigger
    set(gcf, 'color','white','Units', 'Normalized', 'OuterPosition', [0, 0.03, 0.7, 0.8]);
    set(gca, 'FontName', 'Calibri');
    
    plot(sqrt(var_front), ret_front, 'k-', 'LineWidth', 3, 'color', 'b')
hold on;
    plot(sqrt(diag(Q)), mu, 's', 'MarkerSize', 10, 'MarkerFaceColor',[0.9100 0.4100 0.1700]...
        ,'MarkerEdgeColor',[0.9100 0.4100 0.1700])
hold on;
    plot(stdev_random,ret_random, 'o','MarkerSize', 8,'MarkerFaceColor'...
            ,[0 0.5 0], 'MarkerEdgeColor', 'w')
hold off;

% display grid accord to y-axis
    grid on
    ax = gca;
    ax.GridColor = [0.1 0.1 0.1];
    ax.GridLineStyle = '--';
    
% name the graph, x-axis, y-axis and legend
title('EFFICIENT FRONTIER AND DIFFERENT PORTFOLIO STRATEGIES ON MARCH 2015 ','FontSize', 18);
    xlabel('Portfolio Standard Deviation','FontWeight','bold','FontSize',14);
    ylabel('Portfolio Return', 'FontWeight','bold','FontSize',14);
    
lgd = legend('Efficient Frontier', 'Individual Stocks', 'Random Generated Portfolios', ...
    'Location', 'SouthEast');
    lgd.Title.String = 'STRATEGY NAME';
    lgd.FontSize = 14;
            
    