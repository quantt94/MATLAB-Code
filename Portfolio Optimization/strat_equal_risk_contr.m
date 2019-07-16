% EQUAL RISK CONTRIBUTION (ERC PORTFOLIO)
    % Objective: Minimizing deviation of a portfolio w from the equal risk
    % contribution portfolio: min sigma (RC(i) - RC(j))^2 subjected to
    % portfolio weight = 1
    
function [x, cash] = strat_equal_risk_contr(curr_positions, curr_cash, mu, Q, cur_prices)
% number of stocks
n = 20;
curr_prices = cur_prices';

% current value of each component in portfolio (with new prices)
curr_value = curr_positions.*curr_prices;
port_value = sum(curr_value);
curr_w = curr_value./port_value;

% to calculate the target weight of each stock using equal risk
% contributions strategy

global Q A_ineq A_eq


% Equality constraints
A_eq = ones(1,n);
b_eq = 1;

% Inequality constraints
A_ineq = [];
b_ineql = [];
b_inequ = [];
           
% Define initial portfolio weight
w0 = repmat(1.0/n, n, 1);

options.lb = zeros(1,n);       % lower bounds on variables
options.lu = ones (1,n);       % upper bounds on variables
options.cl = [b_eq' b_ineql']; % lower bounds on constraints
options.cu = [b_eq' b_inequ']; % upper bounds on constraints

% Set the IPOPT options
options.ipopt.jac_c_constant        = 'yes';
options.ipopt.hessian_approximation = 'limited-memory';
options.ipopt.mu_strategy           = 'adaptive';
options.ipopt.tol                   = 1e-10;
options.ipopt.print_level           = 0;

% The callback functions
funcs.objective         = @computeObjERC;
funcs.constraints       = @computeConstraints;
funcs.gradient          = @computeGradERC;
funcs.jacobian          = @computeJacobian;
funcs.jacobianstructure = @computeJacobian;
 
%% Run IPOPT
[wsol info] = ipopt(w0',funcs,options);


% Make solution a column vector
if(size(wsol,1)==1)
    w_erc = wsol';
else
    w_erc = wsol;
end

% Compute return, variance and risk contribution for the ERC portfolio
ret_ERC = dot(mu, w_erc)
var_ERC = w_erc'*Q*w_erc;
RC_ERC = (w_erc .* ( Q*w_erc )) / sqrt(w_erc'*Q*w_erc);
 
w_target = w_erc;
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

for i = 1:n;
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
    for i = 1:n;
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
 
