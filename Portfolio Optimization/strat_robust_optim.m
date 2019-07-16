% ROBUST MEAN-VARIANCE OPTIMIZATION
    % Objective: Minimizing variance of portfolio return
    %            Maximizing portfolio expected return
    %            Minimizing Portfolio return estimation error
    %            Constraints: sum of weight = 1 and no short-sales allowed
    
function [x, cash] = strat_robust_optim(curr_positions, curr_cash, mu, Q, cur_prices)
% number of stocks
n = 20;
curr_prices = cur_prices';
r_rf = 0.025;

% current value of each component in portfolio (with new prices)
curr_value = curr_positions.*curr_prices;
port_value = sum(curr_value);
curr_w = curr_value./port_value;

% define initial portfolio: 
% weight = the current weight before balancing accord to the new prices
w0 = curr_w;
ret_init = dot(mu, w0);
var_init = w0' * Q * w0;

% Bounds on variables
lb_rMV = zeros(n,1);
ub_rMV = inf*ones(n,1);

% Required portfolio robustness
var_matr = diag(diag(Q));
% Target portfolio return estimation error is return estimation error of 1/n portfolio
rob_init = w0' * var_matr * w0; % return estimation error of initial portfolio

% target return estimation error is 1 basic point
rob_bnd = 0.01/100; 


% Target portfolio return is risk-free rate assuming for time t = 1/6 year
Portf_Retn = r_rf/6 ;

% Formulate and solve robust mean-variance problem
  
% Objective function
f_rMV  = zeros(n,1);
% Constraints
A_rMV  = sparse([  mu';...
                 ones(1,n)]);
lhs_rMV = [Portf_Retn; 1];
rhs_rMV = [inf; 1];
% Initialize CPLEX environment
cplex_rMV = Cplex('Robust_MV');
% Add objective function and variable bounds
cplex_rMV.addCols(f_rMV, [], lb_rMV, ub_rMV);
% Add constraints
cplex_rMV.addRows(lhs_rMV, A_rMV, rhs_rMV);
% Add quadratic objective
cplex_rMV.Model.Q = 2*Q;
% Add quadratic constraint on return estimation error (robustness constraint)
Qq_rMV = var_matr;
cplex_rMV.addQCs(zeros(size(f_rMV)), Qq_rMV, 'L', rob_bnd, {'qc_robust'});
% Set CPLEX parameters
cplex_rMV.Param.threads.Cur = 4;
cplex_rMV.Param.timelimit.Cur = 60;
cplex_rMV.Param.barrier.qcpconvergetol.Cur = 1e-12; % solution tolerance

% Optimize the problem
cplex_rMV.DisplayFunc = [];
cplex_rMV.solve();   
cplex_rMV.Solution;
    
if(isfield(cplex_rMV.Solution, 'x'))
    w_rMV = cplex_rMV.Solution.x;
    card_rMV = nnz(w_rMV);
    ret_rMV  = dot(mu, w_rMV);
    var_rMV = w_rMV' * Q * w_rMV;
    rob_rMV = w_rMV' * var_matr * w_rMV;
end    
       

% Round near-zero portfolio weights
w_rMV_nonrnd = w_rMV;
w_rMV(find(w_rMV<=1e-6)) = 0;
w_rMV = w_rMV / sum(w_rMV);

w_target = w_rMV;

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

