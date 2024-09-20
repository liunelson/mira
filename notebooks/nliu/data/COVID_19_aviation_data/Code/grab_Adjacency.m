function A = getAstate_multi_update(stateCode,indices,countyIndices,countyAdjacency)
% indices is a vector of fips of counties in the state
nnn=length(stateCode);

n = length(indices);


A = zeros(n,n);


for i=1:n
    target_county=indices(i);
    Whole_ind_county=find(countyIndices==target_county);
    Adjacenty_target_county=countyAdjacency(Whole_ind_county);
    ind_adjacency= find(ismember(indices, Adjacenty_target_county));
    l=length(ind_adjacency);
    for m=1:l
       A(i, ind_adjacency(m))=1;
       A(ind_adjacency(m), i)=1;
    end
end



