% MAXIMUM SHARPE RATIO STRATEGY
    % Objective: Maximizing (mu - Rf)/sigma; subjected to w'i = 1
    % (portfolio weight is 1)
    
function [x, cash] = strat_max_Sharpe(curr_positions, curr_cash, mu, Q, cur_prices)
% number of stocks
N = 20;
curr_prices = cur_prices';
r_rf = 0.025;
% current value of each component in portfolio (with new prices)
curr_value = curr_positions.*curr_prices;
port_value = sum(curr_value);
curr_w = curr_value./port_value;

% to calculate the target weight of each stock using maximum Sharpe ratio
% strategy

% add name to function 
cplex2 = Cplex('max_Sharpe');
cplex2.Model.sense = 'minimize';


% Optimization problem data
lb = zeros(N+1,1);
ub = inf*ones(N+1,1);

% there are 252 trading days within a year
A = (mu - r_rf/252)'; 
b = 1;
c = 0;

% Add objective function and bounds on variables to CPLEX model
cplex2.addCols(zeros(N+1,1), [], lb, ub);

% Add constraints to CPLEX model
cplex2.addRows(b, [A 0], b);
cplex2.addRows(c, [ones(1,N) -1], c);
cplex2.addRows(-inf*ones(N,1), [eye(N) -ones(N,1)], zeros(N,1));

% Add quadratic part of objective function to CPLEX model
cplex2.Model.Q = [2*Q zeros(N,1); zeros(1,N+1)];


% Concurrent algorithm
cplex2.Param.qpmethod.Cur = 6;
% Enable crossover
cplex2.Param.barrier.crossover.Cur = 1;


% Optimize the problem
% cplex2.DisplayFunc = [];
cplex2.solve();

y = cplex2.Solution.x(1:N);
k = cplex2.Solution.x(N+1);

% optimal weight under Max Sharpe Ratio Strategy
w_target = y/k;


% target value of each of the stock in portfolio
target_value = w_target*sum(curr_value);

% the difference between target value and current value = buying and
% selling suggested for each of the stock OR in another way of
% expressing, the total cost of rebalancing the portfolio (without
% transaction costs)
% if negative -->  sell out --> cash inflow
% if positive --> buy more --> cash outflow
adjustment = target_value - curr_value;

% selling/buying units 
x_change = floor(adjustment./curr_prices);

% make sure selling units required would not exceed the amount of stocks we
% have on hands

for i = 1:N;
    if x_change(i) + curr_positions(i) < 0;
        x_change(i) = -1*curr_positions(i);
    end
end

% update new cashflow
cost = x_change.*curr_prices;
transcost = sum(abs(cost).*0.005);
totalcost = sum(cost);

% calculate cash 
cash = curr_cash - totalcost - transcost;
x = curr_positions + x_change;

% if there is not enough cash to rebalance portfolio (negative cash); sell
% part of the portfolio for the total amount equal to the shortage of cash
% without changing the weight of current portfolio

count = 0;
while cash < 0;
    
    % liquidate part of portfolio for cash to rebalance portfolio
    
    % estimated portfolio value after liquidate part of it
    new_port_value = port_value + cash;
    
    % calculate the amount of stocks sold but still keep the current weights
    new_curr_positions = floor((new_port_value.*curr_w)./curr_prices);
    
    % cash on hand recognized by proceed from sale
    proceed = (curr_positions - new_curr_positions).*curr_prices*0.995;
    curr_cash = curr_cash + sum(proceed);
    
    % realized new portfolio value after rounding effect
    port_value = port_value - sum(proceed);
    curr_value = new_curr_positions.*curr_prices;
    
    % rebalance the new portfolio 
    target_value = w_target*port_value;
    adjustment = target_value - curr_value;
    x_change = floor(adjustment./curr_prices);
    for i = 1:N;
    if x_change(i) + new_curr_positions(i) < 0;
        x_change(i) = -1*new_curr_positions(i);
    end
    end
    cost = x_change.*curr_prices;
    transcost = sum(abs(cost).*0.005);
    totalcost = sum(cost);
    cash = curr_cash - totalcost - transcost;
    x = new_curr_positions + x_change;
    
    % counter to observe how many times the portfolio got rebalanced
    count = count + 1;
end
end

