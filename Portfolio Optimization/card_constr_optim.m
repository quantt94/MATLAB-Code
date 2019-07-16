% Compute cardinality-constrained minimum variance portfolio
%

clear all;
clc;

% Add path to CPLEX
addpath('D:/CPLEX/CPLEX1263_x64/cplex/matlab/x64_win64');

% Random data for 10 stocks
n = 10;
Q = randn(n); Q = Q*Q'/1000; % covariance matrix
mu  = rand(1,n)'/100;         % expected return

% Optimization problem data
lb = zeros(n,1);
ub = inf*ones(n,1);
A  = ones(1,n);
b  = 1;

% Define continuous and binary variables
variabtype = [char(ones([1 n])*('C')) char(ones([1 n])*('B'))];

% Compute minimum variance portfolio
cplex = Cplex('min_Variance');
cplex.Model.sense = 'minimize';
cplex.addCols(zeros(2*n,1), [], [lb; zeros(n,1)], [ub; ones(n,1)], variabtype);
cplex.addRows(b, [A zeros(size(A,1),n)], b);
cplex.addRows(n, [zeros(1,n) ones(1,n)], n);
cplex.addRows(-inf*ones(n,1), [eye(n) -eye(n)], zeros(n,1));
cplex.Model.Q = [2*Q zeros(n,n); zeros(n,2*n)];
%cplex.DisplayFunc = []; % disable output to screen

% Compute cardinality-constrained minimum variance portfolios
w_card = [];
for i=1:n
    cplex.Model.lhs(size(A,1)+1) = i;
    cplex.Model.rhs(size(A,1)+1) = i;
    cplex.solve();
    w_card = [w_card cplex.Solution.x(1:n)];
    var_card(i) = w_card(:,i)' * Q * w_card(:,i);
    ret_card(i) = mu' * w_card(:,i);
end

% Plot minimum variance portfolios for cardinality 1, 2, ..., 10
figure(1);
set(gcf, 'color', 'white');
plot([1:n], sqrt(var_card), 'r*', 'MarkerSize', 6)
xlabel('Portfolio cardinality');
ylabel('Portfolio standard deviation');
title('Minimum Variance Portfolios with Cardinality Constraint')