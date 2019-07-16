%% Credit Risk Analysis
% This script follows a simplified version of the JP Morgan Credit Matrics
% approach.  The MathWorks does not endorse or promote the approach and
% places no guarantees on this code; it is simply provided as an example of
% the type of analyses that are possible.
%
% More complete information about CreditMetrics can be found in many
% places, for instance in the technical document located at:
% http://www.defaultrisk.com/pp_model_20.htm

%% 1. Load data
% The needed data is stored in several places in our spreadsheet.

TransitionMatrix  = dataset( ...
    'XLSFile', 'CreditPortfolio.xlsx', ...
    'Sheet', 'Transitions and Rates', ...
    'Range', 'A1:I8', ...
    'ReadObsNames', true, ...
    'ReadVarNames', true);
% This will prove useful:
Ratings = get(TransitionMatrix, 'VarNames');

RecoveryRates = dataset( ...
    'XLSFile', 'CreditPortfolio.xlsx', ...
    'Sheet', 'Transitions and Rates', ...
    'Range', 'A16:C21', ...
    'ReadObsNames', true, ...
    'ReadVarNames', true);
% This will prove useful:
Seniorities = get(RecoveryRates, 'ObsNames');

BondData = dataset( ...
    'XLSFile', 'CreditPortfolio.xlsx', ...
    'Sheet', 'Portfolio Information', ...
    'ReadObsNames', true, ...
    'ReadVarNames', true);
% We can optimize our bond data:
BondData.Rating = ordinal(BondData.Rating, [], Ratings);

BondData.Maturity = datenum(BondData.Maturity, 'mm/dd/yyyy');
BondData.Seniority = ordinal(BondData.Seniority, [], Seniorities);

                        
Rates = dataset( ...
    'XLSFile', 'CreditPortfolio.xlsx', ...
    'Sheet', 'Transitions and Rates', ...
    'Range', 'A10:H14', ...
    'ReadObsNames', true, ...
    'ReadVarNames', true);

CorrelationMatrix = xlsread( ...
    'CreditPortfolio.xlsx', ...
    'Correlations', ...
    'B2:M13');
% We must fill in the upper half of the correlation matrix.
ctrans = CorrelationMatrix';
tempMask = isnan(CorrelationMatrix);
CorrelationMatrix(tempMask) = ctrans(tempMask);
clear ctrans tempMask

CorrelationAlpha = xlsread( ...
    'CreditPortfolio.xlsx', ...
    'Correlations', ...
    'A16');    

numScenarios = xlsread( ...
    'CreditPortfolio.xlsx', ...
    'Risk Analysis', ...
    'A7');

[~,valDate] = xlsread( ...
    'CreditPortfolio.xlsx', ...
    'Risk Analysis', ...
    'B7');
valDate = datenum(valDate, 'mm/dd/yyyy');

PercentForVaR = xlsread( ...
    'CreditPortfolio.xlsx', ...
    'Risk Analysis', ...
    'A2:A4');
PercentForVaR = 100 * PercentForVaR;

numAssets = size(BondData,1);

%% 2. Compute original value of portfolio
dates = [valDate; valDate+365; valDate+2*365; valDate+3*365];

% The |prbyzero| function in the Financial Toolbox allows us to value these
% bonds.

for idx = 1 : length(Ratings)-1
    ratingsMask = BondData.Rating == Ratings{idx};
    BondData.ValuePerBond(ratingsMask) = prbyzero( ...
        [BondData.Maturity(ratingsMask) BondData.Coupon(ratingsMask)],...
        valDate, Rates.(Ratings{idx}), dates);
end

OriginalPortfolioValue = BondData.ValuePerBond' * BondData.NumberOfBonds;

%% 3. Compute ratings thresholds
% We now need to calculate the ratings thresholds. These represent the
% number of standard deviations an asset price has to move before the
% parent company suffers a rating change.
% 
% Using the Asset Value Model, we will compute for each of the ratings the
% Zdefault, Zccc, Zb, Zbb, Zbbb, Za, Zaa, and Zaaa thresholds. We will
% start with Zdefault:
%
% $$ P_{Default}=P(R<Z_{Default})=\phi(Z_{Default}\mid\sigma) $$
%
% and then iteratively compute the other thresholds using the following
% formula :
%
% $$ P_{CCC}=P(Z_{Default}<R<Z_{CCC})=\phi(Z_{CCC}\mid\sigma)-\phi(Z_{Default} \mid\sigma) $$

% Cumulatively add the transition probabilities, starting at the lower end.
ThresholdMatrix = fliplr(cumsum(fliplr(double(TransitionMatrix)), 2));
% Turn these into thresholds.
ThresholdMatrix = norminv(ThresholdMatrix);
% Convert into a dataset array with a similar format as |TransitionMatrix|.
ThresholdMatrix = replacedata(TransitionMatrix, ThresholdMatrix); 

%% 4. Generate scenarios
% We will generate a large number of scenarios. To do this, we create a set
% of random numbers correlated by industry within each of the scenarios.
% For each bond we will then add an idiosyncratic random number
% corresponding to the variation within each industry.  The parameter
% governing the relative importance of these two random processes is
% |CorrelationAlpha|, which for the sake of simplicity we have assumed to
% be constant for each bond in our portfolio.

IndustryRandom = mvnrnd( ...
    zeros(size(CorrelationMatrix, 1), 1), ...
    CorrelationMatrix, ...
    numScenarios)';

% By adjusting the weights properly in terms of |CorrelationAlpha|, we
% ensure that the marginal distribution of each company's performance is a
% standard, normally-distributed number.  This mathes the requirements of
% our ratings thresholds.
BondRandom = CorrelationAlpha*IndustryRandom(BondData.Industry,:) + ...
    sqrt(1-CorrelationAlpha.^2)*randn(numAssets,numScenarios);

%% 5. Find rating of each bond at the end of the simulations
% We use the threshold levels and the generated scenarios. For each
% scenario, find how the companies fared, and then give them a new rating
% if they have crossed a threshold.
NewRatings = ones(numAssets,numScenarios);

for i = 1 : numAssets
    FirmRating = BondData.Rating(i);
    Thresholds = [double(ThresholdMatrix(char(FirmRating),:)) -Inf];
    for j=2:numel(Thresholds)
        RatingIndex = find(BondRandom(i,:) > Thresholds(j) & BondRandom(i,:) < Thresholds(j-1));
        NewRatings(i,RatingIndex) = j-1;
    end
end

% We no longer need this large variable.
clear BondRandom

NewRatings = ordinal(NewRatings, Ratings);

%% 6. Value bonds based on the simulated ratings
% Now we need to value the bonds for each of the scenarios. This is similar
% to step 2: we price each bond if the bond has not defaulted, and if it
% has, then we calculate how much we expect to recover.
%
% It will be a significant time saver if we pre-calculate the non-default
% prices.  This will save us from having to re-calculate the prices
% multiple times.
NondefaultPrices = prbyzero([BondData.Maturity BondData.Coupon], ...
    valDate, double(Rates), dates);

NewPrices = zeros(numAssets,numScenarios);

for Bond = 1 : numAssets
    for idx = 1 : length(Ratings)
        
        % Consider nondefault cases:
        if idx <= 7; 
            ratingsMask = NewRatings(Bond,:) == Ratings{idx};
            NewPrices(Bond,ratingsMask) = NondefaultPrices(Bond,idx);
        
        % Consider default cases:
        else
            Seniority = char(BondData{Bond,'Seniority'});
            DistribParameters = double(RecoveryRates(Seniority,{'Mean','STD'}));
            defaultMask = NewRatings(Bond,:) == 'Default';
            NbToDraw = sum(defaultMask);
            NewPrices(Bond, defaultMask) = ...
                100 * betarndms(DistribParameters(1),DistribParameters(2),...
                1,NbToDraw);
        end
    end
end

%% 7. Compute values of simulated portfolios
% The |bsxfun| function is designed to efficiently expand a vector (like
% the 1311x1 "number of bonds" so that it can be multiplied
% element-by-element against a matrix (like the 1311x10000 "New Prices").
NewPrices = bsxfun(@times, NewPrices, BondData.NumberOfBonds);
NewPortfolioValues = sum(NewPrices,1);

%% 8. Display the distribution of possible values for the portfolio

Losses = OriginalPortfolioValue - NewPortfolioValues;

VaR = prctile(Losses, PercentForVaR);
relativeVaR = VaR / OriginalPortfolioValue;

CVaR = zeros(size(VaR));

for i = 1:length(PercentForVaR)
    CVaR(i) = mean(Losses(Losses >= VaR(i)));
end

relativeCVaR = CVaR / OriginalPortfolioValue;

figure, hold on,
set(gcf,'Position',[100 100 800 350])
hist(Losses,100);
l1 = plot([VaR(2) VaR(2)], ylim, 'b', 'LineWidth', 2);
l2 = plot([CVaR(2) CVaR(2)], ylim, 'r', 'LineWidth', 2);
set(gca, 'FontSize', 16, 'FontName', 'Calibri', ...
    'YTickLabel','')
xlabel('Losses (in USD)')
legend([l1 l2], ...
    [int2str(PercentForVaR(2)) '% VaR: ' cur2str(VaR(2),0)], ...
    [int2str(PercentForVaR(2)) '% CVaR: ' cur2str(CVaR(2),0)])
hold off