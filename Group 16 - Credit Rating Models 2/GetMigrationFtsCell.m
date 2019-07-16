function migrFtsCell = GetMigrationFtsCell(id,date,numRating)

nCo = 0; % number of companies in the data set, to be determined here
currID = id(1)-1; % aux array to keep track of IDs
nRecords = length(id); % number of records in the data
startIdx = zeros(nRecords+1,1); % Array of "pointers": row in the data set
                                % where each company starts
for i = 1:nRecords
   if (id(i)~=currID)
      nCo = nCo + 1;
      startIdx(nCo) = i;
      currID = id(i);
   end
end
startIdx(nCo+1) = nRecords + 1; % This last entry is only a trick to
                                % simplify the code below
migrFtsCell = cell(nCo,1);
for k = 1:nCo
   migrFtsCell{k} = fints(date(startIdx(k):startIdx(k+1)-1),...
      numRating(startIdx(k):startIdx(k+1)-1),...
      'rating', 0,...
      strcat('Rating migration history, co. ',num2str(id(startIdx(k)))));
end

end