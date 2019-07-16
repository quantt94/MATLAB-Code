% EQUALLY WEIGHTED STRATEGY 
    % Objective: asset weights are selected as w(i) = 1/n where n = 20 in
    % this case (total number of stocks in portfolio)
    
function [x, cash] = strat_equally_weighted(curr_positions, curr_cash, mu, Q, cur_prices);

N = 20;
curr_prices = cur_prices';

% current value of each component in portfolio (with new prices)
curr_value = curr_positions.*curr_prices;
port_value = sum(curr_value);
curr_w = curr_value./port_value;

% target value of each component under equally weighted strategy
target_value = (sum(curr_value)/N)*ones(N,1);

% the difference between target value and current value = buying and
% selling amount suggested for each of the stock OR in another way of
% expressing, the total cost of rebalancing the portfolio (without
% transaction costs)
% if negative -->  sell out --> cash inflow
% if positive --> buy more --> cash outflow
adjustment = target_value - curr_value;

% selling/buying units 
x_change = floor(adjustment./curr_prices);

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
    
    % realized new portfolio value after rounding effect; sum(curr_value)
    % should equal to port_value 
    port_value = port_value - sum(proceed);
    curr_value = new_curr_positions.*curr_prices;
    
    % rebalance the new portfolio 
    target_value = (port_value/N)*ones(N,1);
    adjustment = target_value - curr_value;
    x_change = floor(adjustment./curr_prices);
    cost = x_change.*curr_prices;
    transcost = sum(abs(cost).*0.005);
    totalcost = sum(cost);
    cash = curr_cash - totalcost - transcost;
    x = new_curr_positions + x_change;
    
    % counter to observe how many times the portfolio got rebalanced
    count = count + 1;
    
end
end
    
    
    
    



