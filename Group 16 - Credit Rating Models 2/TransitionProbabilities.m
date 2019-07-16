%% Estimating Transition Probabilities
clc;
clear;
%% 1. Import the data
% MATLAB's dataset arrays (available within Statistics Toolbox) provide a
% convenient way to load in this data, which is a combination of numeric
% (ID numbers) and textual (dates and ratings) information. 
rawData = dataset('File','MigrationHistory.txt','VarNames',...
                  {'ID','Date','Rating'},'ReadVarNames',false,...
                  'Delimiter',',');

% We can examine a piece of this dataset with commands like the following:
% firstCompany = rawData(rawData.ID == 10283,:) 

%% 2. Convert raw data using ordinal arrays
% Ordinal arrays offer several benefits for storing ordered, label-like
% information (like bond ratings).  For example, they occupy a much smaller
% memory footprint, and operations on ordinal array will execute faster
% than those on cell arrays.
rawData.Rating = ordinal(rawData.Rating,[],...
                         {'AAA','AA','A','BBB','BB','B','CCC','D'});
                     
%% FIRST CHALLENGE: DATA MANAGEMENT

%% 3. Create financial time series for first company
% Financial time series objects are particularly useful for time and data
% manupulations.  In our case, we need to resample our rating data to
% one-year intervals.  Let us look at a single company as an example.
myFts = fints(rawData.Date(1:4), double(rawData.Rating(1:4)),...
              'rating', 0);

%% 4. Define time window

startDate = '12/31/1986';
endDate   = '12/31/2000';
dates     = cfdates(datenum(startDate)-1,datenum(endDate),1)';

%% 5. Fill missing values in time series
% The |fillts| function fills in the resampled annual observations.
myFts = fillts(myFts,'zero',dates, 1);

%% 6. Extract relevant dates from time series

myFts = myFts(datestr(dates));

%% SECOND CHALLENGE: PERFORMANCE

%% 7. Prepare for computing migration matrix
% Let us repeat the above steps for all of the historical companies.
startDate   = '12/31/1986';
endDate     = '12/31/2000';
dates       = cfdates(datenum(startDate)-1,datenum(endDate),1)';
migrFtsCell = GetMigrationFtsCell(rawData.ID,rawData.Date,...
                                  double(rawData.Rating));
                                  
%% 8. Compute migration matrix
migrMat = zeros(size(migrFtsCell,1),length(dates));

% MATLAB's |parfor|, or parallel FOR-loop, allows us to easily take
% advantage of multiple processors if they are available.  In a
% single-processor environment, this will run as a simple FOR-loop, but in
% a parallel environment, MATLAB will automatically send out iterations of
% the loop to workers as they become available.  This results in faster
% computing.
tic
parfor k = 1:size(migrFtsCell,1)
   ftsDummy = fillts(migrFtsCell{k},'zero',dates,1);
   [~,idx,~] = intersect(ftsDummy.dates,dates);
   migrMat(k,:) = fts2mat(ftsDummy(idx));
end
toc

%% 9. Compute transition probability matrix

nRatings  = 8;
transProb = GetTransProb(migrMat,nRatings);

%% VISUALIZATION - SOCIAL COMPUTING
% We can use a visualization function from MATLAB Central to visualize the
% transition matrix as a heatmap.

ratingsList = {'AAA'; 'AA'; 'A'; 'BBB'; 'BB'; 'B'; 'CCC'; 'D'};
heatmap(transProb*100, ratingsList, ratingsList, '%0.2f%%',...
        'Colormap','red','UseLogColorMap',true,'Colorbar',true);
