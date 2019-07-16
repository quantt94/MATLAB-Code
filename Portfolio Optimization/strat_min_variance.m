% MINIMUM VARIANCE STRATEGY 
    % Objective: Minimizing porfolio variance - w'Qw subjected to w'i = 1
    % (portfolio weight = 1)
    
function [x, cash] = strat_min_variance(curr_positions, curr_cash, mu, Q, cur_prices)
% number of stocks
N = 20;
curr_prices = cur_prices';
% current value of each component in portfolio (with new prices)
curr_value = curr_positions.*curr_prices;
port_value = sum(curr_value);
curr_w = curr_value./port_value;

% to calculate target weight of each stock by 
% minimizing variance 

lb = zeros(N,1);
ub = inf*ones(N,1);
A  = ones(1,N);
b  = 1;

% add name to function
cplex1 = Cplex('min_Variance');
% Add objective function and bounds on variables to CPLEX model
cplex1.addCols(zeros(N,1), [], lb, ub); 
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

% Read results from minimum variance portfolio
w_target = cplex1.Solution.x;
var_minVar = w_target' * Q * w_target;
ret_minVar = mu' * w_target;

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

