function transProb = GetTransProb(migrMat,nRatings)

nInit = zeros(nRatings,1); % companies that start a period as rating i
nTrans = zeros(nRatings,nRatings); % companies that migrate from i to j

[nCo,T] = size(migrMat);

for k = 1:nCo
    for j = 1:T-1
        nInit(migrMat(k,j)) = nInit(migrMat(k,j)) + 1;
        nTrans(migrMat(k,j),migrMat(k,j+1)) =...
            nTrans(migrMat(k,j),migrMat(k,j+1)) + 1;
    end
end

transProb = zeros(size(nTrans));
nonEmpty = find(nInit>0); % Avoid dividing by zero for unvisited ratings
transProb(nonEmpty,:) = diag(nInit(nonEmpty)) \ nTrans(nonEmpty,:);

end