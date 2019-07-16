clc;
clear all;
format long;

input_file_prices  = 'Daily_closing_prices.csv';

% read data
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


% compute mu and Q
hist_return = data_prices(2:end,:) ./ data_prices(1:end-1,:) - 1;
mu = mean(hist_return)';
Q = cov(hist_return);

n = 20;
lb = zeros(n,1);
ub = inf*ones(n,1);
A = ones(1,n);
b = 1;




m=50;
[M N]   = size(hist_return);
w_minVar = zeros(20, m);
var_minVar = zeros(1, m);
ret_minVar = zeros(1, m);
w_maxRet = zeros(20, m);
var_maxRet = zeros(1, m);
ret_maxRet = zeros(1, m);
for i=1:m
    % Simulation
    simRets = mvnrnd(mu, Q, M);
    %simRets(:,18) = simRets(:,18)*1.001;
    muSim{i} = mean(simRets)';
    covSim{i} = cov(simRets);
    
    % compute minimum variance portfolio
    cplex1 = Cplex('min_Variance');
    cplex1.addCols(zeros(n,1), [], lb, ub);
    cplex1.addRows(b, A, b);
    cplex1.Model.Q = 2*covSim{i};
    cplex1.Param.qpmethod.Cur = 6;
    cplex1.Param.barrier.crossover.Cur = 1;
    cplex1.DisplayFunc = [];
    cplex1.solve();
    w_minVar(:,i) = cplex1.Solution.x;
    var_minVar(i) = w_minVar(:,i)' * covSim{i} * w_minVar(:,i);
    ret_minVar(i) = muSim{i}'* w_minVar(:,i);
    
    % compute maximum return portfolio
    cplex2 = Cplex('max_Return');
    cplex2.Model.sense = 'maximize';
    cplex2.addCols(muSim{i}, [], lb, ub);
    cplex2.addRows(b, A, b);
    cplex2.Param.qpmethod.Cur = 6;
    cplex2.Param.barrier.crossover.Cur = 1;
    cplex2.DisplayFunc = [];
    cplex2.solve();
    w_maxRet(:,i) = cplex2.Solution.x;
    var_maxRet(i) = w_maxRet(:,i)' * covSim{i} * w_maxRet(:,i);
    ret_maxRet(i) = muSim{i}'* w_maxRet(:,i);

    %compute efficient frontier
    targetRet = linspace(ret_minVar(i), ret_maxRet(i), 30);
    
    cplex3 = Cplex('Efficient_Frontier');
    cplex3.addCols(zeros(n,1), [], lb, ub);
    cplex3.addRows(targetRet(1), muSim{i}', inf);
    cplex3.addRows(b, A, b);
    cplex3.Model.Q = 2*covSim{i};
    cplex3.Param.qpmethod.Cur = 6;
    cplex3.Param.barrier.crossover.Cur = 1;
    cplex3.DesplayFunc = [];
    
    w_front = [];
    for j=1:length(targetRet)
        cplex3.Model.lhs(1) =  targetRet(j);
        cplex3.solve();
        w_front = [w_front cplex3.Solution.x];
        var_front(j) = w_front(:,j)' * covSim{i} * w_front(:,j);
        ret_front(j) = muSim{i}' * w_front(:,j);
    end
    
    sim_var_front(i,:) = var_front;
    sim_ret_front(i,:) = ret_front;
    MarkoWR{i} = w_front;
    
end

figure(1)
subplot(1,2,1);
for k = 1:m
    color = rand(1,3);
    plot(sqrt(sim_var_front(k,:)),sim_ret_front(k,:),'Color',color)
    hold on
    plot(sqrt(diag(covSim{k})),muSim{k},'r.','Color',color)
end
set(gca,'Box','on','LineWidth', 1.5, 'FontSize',12');
title('Resampled frontiers');
ylabel('returns');
xlabel('volatility');
grid on



% Compute minimum variance portfolio
cplex1 = Cplex('min_Variance');
cplex1.addCols(zeros(n,1), [], lb, ub);
cplex1.addRows(b, A, b);
cplex1.Model.Q = 2*Q;
cplex1.Param.qpmethod.Cur = 6;
cplex1.Param.barrier.crossover.Cur = 1;
cplex1.DisplayFunc = [];
cplex1.solve();

% Display minimum variance portfolio
w_minVar = cplex1.Solution.x;
var_minVar = w_minVar' * Q * w_minVar;
ret_minVar = mu'* w_minVar;


% Compute maximum retrun portfolio
cplex2 = Cplex('max_Return');
cplex2.Model.sense = 'maximize';
cplex2.addCols(mu, [], lb, ub);
cplex2.addRows(b, A, b);
%cplex2.Model.Q = 2*Q;
cplex2.Param.qpmethod.Cur = 6;
cplex2.Param.barrier.crossover.Cur = 1;
cplex2.DisplayFunc = [];
cplex2.solve();

% Display maximum return portfolio
w_maxRet = cplex2.Solution.x;
var_maxRet = w_maxRet' * Q * w_maxRet;
ret_maxRet = mu'* w_maxRet;


% Target returns
targetRet = linspace(ret_minVar, ret_maxRet, 30);

% Compute efficient frontier
cplex3 = Cplex('Efficient_Frontier');
cplex3.addCols(zeros(n,1), [], lb, ub);
cplex3.addRows(targetRet(1), mu', inf);
cplex3.addRows(b, A, b);
cplex3.Model.Q = 2*Q;
cplex3.Param.qpmethod.Cur = 6;
cplex3.Param.barrier.crossover.Cur = 1;
cplex3.DesplayFunc = [];

w_front = [];
for i=1:length(targetRet)
    cplex3.Model.lhs(1) =  targetRet(i);
    cplex3.solve();
    w_front = [w_front cplex3.Solution.x];
    var_front(i) = w_front(:,i)' * Q * w_front(:,i);
    ret_front(i) = mu' * w_front(:,i);
end
subplot(1,2,2);
set(gcf, 'color', 'white');
plot(sqrt(var_front), ret_front, 'k-', 'LineWidth', 3)
hold on;
plot(sqrt(var_minVar), ret_minVar, 'rd', 'MarkerSize', 6)
hold on;
plot(sqrt(var_maxRet), ret_maxRet, 'ms', 'MarkerSize', 6)
hold on;
plot(sqrt(diag(Q)), mu, 'b.', 'MarkerSize', 18)
xlabel('Standard deviation');
ylabel('Expected return');
title('Efficient Frontier');
legend('efficient frontier', 'minimum variance portfolio', 'maximum return portfolio', 'individual stocks', 'Location', 'SouthWest');

% calculate average
aveW = zeros(20, 30);
for t = 1:30
    for  k = 1:m
        aveW(:,t) =  aveW(:,t) + MarkoWR{k}(:,t);
    end
end
aveW = 1/m * aveW;
for s = 1:30
    var_ave(s) = aveW(:,s)' * Q * aveW(:,s);
    ret_ave(s) = mu' * aveW(:,s);
end
figure(2)
set(gcf, 'color', 'white');
plot(sqrt(var_ave), ret_ave, 'b', 'LineWidth', 3)
hold on;
plot(sqrt(var_front), ret_front, 'r--', 'LineWidth', 3)
xlabel('Standard deviation');
ylabel('Expected return');
title('Average Efficient Frontier');

figure(3)
subplot(1,2,1)
area(w_front')
title('Markowitz Portfolio Weights');
set(get(gcf,'Children'),'YLim',[0 1]);
xlabel('portfolio')
ylabel('asset weight')
subplot(1,2,2)
area(aveW')
title('Resampled Portfolio Weights');
set(get(gcf,'Children'),'YLim',[0 1]);
xlabel('portfolio');
ylabel('asset weight');
