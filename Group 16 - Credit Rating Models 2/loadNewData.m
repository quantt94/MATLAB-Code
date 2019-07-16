function [newdata, IDs] = loadNewData(fileToRead1)
%LOADNEWDATA(FILETOREAD1)
%  Imports data from the specified file
%  FILETOREAD1:  file to read

%  Auto-generated by MATLAB on 05-Feb-2010 10:43:27

% Import the file
sheetName = 'Credit Rating';
[numbers, headers] = xlsread(fileToRead1, sheetName);


% Filter this data to only include the significant predictors.  Note that
% they may be in a different order than for our training data set, so we
% will use the header information to find the proper columns:
reInd   = find(strcmp(headers, 'RE_TA'));
mveInd  = find(strcmp(headers, 'MVE_BVTD'));
indInd  = find(strcmp(headers, 'Industry'));

newdata = numbers(:, [reInd mveInd indInd]);

% We will also read in the company IDs, which may be useful for reporting
% purposes.
idInd   = strcmp(headers, 'ID');
IDs     = numbers(:, idInd);
