function [LCS,tcoeff] = CreateLCS(p,pc,index,K)

count = pc.Count;

KNN = findNearestNeighbors(pc,p.Location,K);
    
points = zeros(numel(KNN),3);

i = 1:numel(KNN);
p_i = select(pc,KNN(i));
points(i,:) = p_i.Location;

npoints = zscore(points);

coeff = pca(npoints);

tcoeff = transpose(coeff);

%biplot(coeff);

LCS(1,:) = [ 0 0 0 ];

j = 1:count;

pj = select(pc,j(j<index));

temp = transpose(pj.Location - p.Location);

LCS(j(j<index)+1,:) = transpose(tcoeff*temp);

pj = select(pc,j(j>index));

temp = transpose(pj.Location - p.Location);

LCS(j(j>index),:) = transpose(tcoeff*temp);

end