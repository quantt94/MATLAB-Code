clc;
clear all;
format long

warning('off','all')
warning

% Input files
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

% Find the number of trading days in Nov-Dec 2014 and
% compute expected return and covariance matrix for period 1
day_ind_start0 = 1;
day_ind_end0 = length(find(dates_array(:,1)==2014));
cur_returns0 = data_prices(day_ind_start0+1:day_ind_end0,:) ./ data_prices(day_ind_start0:day_ind_end0-1,:) - 1;
mu = mean(cur_returns0)';
Q = cov(cur_returns0);

% Remove datapoints for year 2014
data_prices = data_prices(day_ind_end0+1:end,:);
dates_array = dates_array(day_ind_end0+1:end,:);
dates = dates(day_ind_end0+1:end,:);

% Initial positions in the portfolio
init_positions = [5000 950 2000 0 0 0 0 2000 3000 1500 0 0 0 0 0 0 1001 0 0 0]';

% Initial value of the portfolio
init_value = data_prices(1,:) * init_positions;
fprintf('\nInitial portfolio value = $ %10.2f\n\n', init_value);

% Initial portfolio weights
w_init = (data_prices(1,:) .* init_positions')' / init_value;

% Number of periods, assets, trading days
N_periods = 6*length(unique(dates_array(:,1))); % 6 periods per year
N = length(tickers);
N_days = length(dates);

% Annual risk-free rate for years 2015-2016 is 2.5%
r_rf = 0.025;
% Annual risk-free rate for years 2008-2009 is 4.5%
r_rf2008_2009 = 0.045;

% Number of strategies
strategy_functions = {'strat_buy_and_hold' 'strat_equally_weighted' 'strat_min_variance' 'strat_max_Sharpe' 'strat_equal_risk_contr' 'strat_lever_equal_risk_contr' 'strat_robust_optim'};
strategy_names     = {'Buy and Hold' 'Equally Weighted Portfolio' 'Minimum Variance Portfolio' 'Maximum Sharpe Ratio Portfolio' 'Equal Risk Contributions Portfolio' 'Leveraged Equal Risk Contributions Portfolio' 'Robust Optimization Portfolio'};
N_strat = 7; % comment this in your code
%N_strat = length(strategy_functions); % uncomment this in your code
fh_array = cellfun(@str2func, strategy_functions, 'UniformOutput', false);

for (period = 1:N_periods)
   % Compute current year and month, first and last day of the period
   if(dates_array(1,1)==15)
       cur_year  = 15 + floor(period/7);
   else
       cur_year  = 2015 + floor(period/7);
   end
   cur_month = 2*rem(period-1,6) + 1;
   day_ind_start = find(dates_array(:,1)==cur_year & dates_array(:,2)==cur_month, 1, 'first');
   day_ind_end = find(dates_array(:,1)==cur_year & dates_array(:,2)==(cur_month+1), 1, 'last');
   fprintf('\nPeriod %d: start date %s, end date %s\n', period, char(dates(day_ind_start)), char(dates(day_ind_end)));

   % Prices for the current day
   cur_prices = data_prices(day_ind_start,:);

   % Execute portfolio selection strategies
   for(strategy = 1:N_strat)

      % Get current portfolio positions
      if(period==1)
         curr_positions = init_positions;
         curr_cash = 0;
         portf_value{strategy} = zeros(N_days,1);
      else
         curr_positions = x{strategy,period-1};
         curr_cash = cash{strategy,period-1};
      end

      % Compute strategy
      [x{strategy,period} cash{strategy,period}] = fh_array{strategy}(curr_positions, curr_cash, mu, Q, cur_prices);

      % Verify that strategy is feasible (you have enough budget to re-balance portfolio)
      % Check that cash account is >= 0
      % Check that we can buy new portfolio subject to transaction costs

      %%%%%%%%%%% Insert your code here %%%%%%%%%%%%
      % to calculate the allocation of each of the stock (the weight
      % corresponding to total portfolio value for each of the period)
      allocation{strategy,period} = (cur_prices'.*x{strategy,period})./(cur_prices*x{strategy,period});


      % Compute portfolio value
      portf_value{strategy}(day_ind_start:day_ind_end) = data_prices(day_ind_start:day_ind_end,:) * x{strategy,period} + cash{strategy,period};

      fprintf('   Strategy "%s", value begin = $ %10.2f, value end = $ %10.2f\n', char(strategy_names{strategy}), portf_value{strategy}(day_ind_start), portf_value{strategy}(day_ind_end));

   end
      
   % Compute expected returns and covariances for the next period
   cur_returns = data_prices(day_ind_start+1:day_ind_end,:) ./ data_prices(day_ind_start:day_ind_end-1,:) - 1;
   mu = mean(cur_returns)';
   Q = cov(cur_returns);
   
end

% Plot results
% figure(1);
%%%%%%%%%%% Insert your code here %%%%%%%%%%%%

% PLOT DAILY PORTFOLIO VALUE FOR DIFFERENT STRATEGIES 

% Convert dates to string for the use in time series plot
date = datestr(dates);

% Create a loop to plot multi-line within a graph 
% Set colors to multiple lines
colvec = [0 0.5 1;
          1 0.2 0.2;
          0 0.6 0.3;
          1 0.5 0;
          0.2 0.8 0.8;
          0 0 1;
          0.2 0 0];
for i = [1:5 7];
    ts(i) = timeseries(portf_value{1,i},date);
    ts(i).TimeInfo.Format = 'mmm dd, yy';
    plot(ts(i),'Color', colvec(i,:),'LineWidth', 1.2);
    hold on
    
end
    % set the limit on x axis
    tstart = datetime(2015, 01, 01);
    tend = datetime(2016, 12, 30);
    xlim([tstart tend]);
    
    % set the plotbox aspect ratio 
    pbaspect([1 0.514 0.514]);
    
    % display grid accord to y-axis
    grid on
    ax = gca;
    ax.GridColor = [0.1 0.1 0.1];
    ax.GridLineStyle = '--';
    
    % set the figure size to be bigger
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.03, 0.7, 0.8]);
    set(gca, 'FontName', 'Calibri');
    
    % name the graph, x-axis, y-axis and legend
    title('PORTFOLIO VALUE OVER THE YEARS 2015 AND 2016','FontSize', 18);
    xlabel('Dates','FontWeight','bold','FontSize',14);
    ylabel('Portfolio Value (in Millions)', 'FontWeight','bold','FontSize',14);
    lgd = legend('Buy and Hold','Equally Weighted','Minimum Variance','Maximum Sharpe Ratio','Equal Risk Contributions','Robust Mean-Variance Optimization','location','northwest');
    lgd.Title.String = 'STRATEGY NAME';
    lgd.FontSize = 14;

    % plot portfolio value for LEVERAGED EQUAL RISK CONTRIBUTIONS
    figure;
    ts6 = timeseries(portf_value{1,6},date);
    ts6.TimeInfo.Format = 'mmm dd, yy';
    plot(ts6,'Color', colvec(2,:),'LineWidth', 1.2);
    % set the limit on x axis
    tstart = datetime(2015, 01, 01);
    tend = datetime(2016, 12, 30);
    xlim([tstart tend]);
    
    % set the plotbox aspect ratio 
    pbaspect([1 0.514 0.514]);
    
    % display grid accord to y-axis
    grid on
    ax = gca;
    ax.GridColor = [0.1 0.1 0.1];
    ax.GridLineStyle = '--';
    
    % set the figure size to be bigger
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.03, 0.7, 0.8]);
    set(gca, 'FontName', 'Calibri');
    
    % name the graph, x-axis, y-axis and legend
    title('PORTFOLIO VALUE OVER THE YEARS 2015 AND 2016','FontSize', 18);
    xlabel('Dates','FontWeight','bold','FontSize',14);
    ylabel('Portfolio Value (in Millions)', 'FontWeight','bold','FontSize',14);
    lgd = legend('Leveraged Equal Risk Contributions','location','northwest');
    lgd.Title.String = 'STRATEGY NAME';
    lgd.FontSize = 14;
    
    %% Plot chart for dynamic change in position 

% Create vector of securities' names
for i = 2:21;
    name{i-1} = headers{1,1}{i,1};
end
name = name';
name = string(name);

% MINIMUM VARIANCE STRATEGY

% Combine the allocation changes into one matrix
    for i = 1:12;
        minvar(:,i) = allocation{3,i};
    end    
               
% Create multiple lines 
figure;
cmap = jet(20);
    for i = 1:20;
        plot(minvar(i,:), '-s','Color', cmap(i,:),'LineWidth',2.0);
        legendInfo{i} = name(i);
        hold on;
    end
    % set x-axis and y-axis limit
    xlim([1 12]);
    ylim([-0.1 1]);
    
    % set the plotbox aspect ratio 
    pbaspect([1 0.514 0.514]);
    
    % display grid accord to y-axis
    grid on
    ax = gca;
    ax.GridColor = [0.1 0.1 0.1];
    ax.GridLineStyle = '--';
    
    % set the figure size to be bigger
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.03, 0.7, 0.8]);
    set(gca, 'FontName', 'Calibri','Projection','perspective');
    
    % name the graph, x-axis, y-axis and legend
    title('DYNAMIC CHANGES IN PORTFOLIO ALLOCATIONS (Minimum Variance Strategy)','FontSize', 18);
    xlabel('Periods (bi-monthly)','FontWeight','bold','FontSize',14);
    ylabel('Allocation of each stock in Portfolio', 'FontWeight','bold','FontSize',14);
    lgd = legend(legendInfo,'location','westoutside');
    lgd.Title.String = 'STOCK NAME';
    lgd.FontSize = 14;
    
% MAXIMUM SHARPE RATIO STRATEGY

% Combine the allocation changes into one matrix
    for i = 1:12;
        maxsharpe(:,i) = allocation{4,i};
    end    
               
% Create multiple lines 
figure;
cmap = jet(20);
    for i = 1:20;
        plot(maxsharpe(i,:), '-*','Color', cmap(i,:),'LineWidth',2.0);
        legendInfo{i} = name(i);
        hold on;
    end
    % set x-axis and y-axis limit
    xlim([1 12]);
    ylim([-0.1 1]);
    
    % set the plotbox aspect ratio 
    pbaspect([1 0.514 0.514]);
    
    % display grid accord to y-axis
    grid on
    ax = gca;
    ax.GridColor = [0.1 0.1 0.1];
    ax.GridLineStyle = '--';
    
    % set the figure size to be bigger
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.03, 0.7, 0.8]);
    set(gca, 'FontName', 'Calibri','Projection','perspective');
    
    % name the graph, x-axis, y-axis and legend
    title('DYNAMIC CHANGES IN PORTFOLIO ALLOCATIONS (Maximum Sharpe Ratio Strategy)','FontSize', 18);
    xlabel('Periods (bi-monthly)','FontWeight','bold','FontSize',14);
    ylabel('Allocation of each stock in Portfolio', 'FontWeight','bold','FontSize',14);
    lgd = legend(legendInfo,'location','westoutside');
    lgd.Title.String = 'STOCK NAME';
    lgd.FontSize = 14;   
    
% ROBUST MEAN-VARIANCE OPTIMIZATION PORTFOLIO

% Combine the allocation changes into one matrix
    for i = 1:12;
        maxsharpe(:,i) = allocation{7,i};
    end    
               
% Create multiple lines 
figure;
cmap = jet(20);
    for i = 1:20;
        plot(maxsharpe(i,:), '-*','Color', cmap(i,:),'LineWidth',2.0);
        legendInfo{i} = name(i);
        hold on;
    end
    % set x-axis and y-axis limit
    xlim([1 12]);
    ylim([-0.1 1]);
    
    % set the plotbox aspect ratio 
    pbaspect([1 0.514 0.514]);
    
    % display grid accord to y-axis
    grid on
    ax = gca;
    ax.GridColor = [0.1 0.1 0.1];
    ax.GridLineStyle = '--';
    
    % set the figure size to be bigger
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.03, 0.7, 0.8]);
    set(gca, 'FontName', 'Calibri','Projection','perspective');
    
    % name the graph, x-axis, y-axis and legend
    title('DYNAMIC CHANGES IN PORTFOLIO ALLOCATIONS (Robust Mean-Variance Optimization)','FontSize', 18);
    xlabel('Periods (bi-monthly)','FontWeight','bold','FontSize',14);
    ylabel('Allocation of each stock in Portfolio', 'FontWeight','bold','FontSize',14);
    lgd = legend(legendInfo,'location','westoutside');
    lgd.Title.String = 'STOCK NAME';
    lgd.FontSize = 14;
    
