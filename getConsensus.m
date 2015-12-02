function [ consensus ] = getConsensus( cloud1 )
%getting consensus sequence out of genetic cloud


nucs = {'A', 'C', 'G', 'T', 'N', '-'};
%to initialise consensus, error if still there at the end
consensus(size(cloud1,2)) = 'E';
for i=1:size(cloud1,2)
    [~,pos] =  max( cloud1(:,i) );
    consensus(i) = nucs{pos} ;
    %if sum(cloud1(:,i)) < mean(sum(cloud1))/4
    if sum(cloud1(:,i)) < 10
        consensus(i) = 'N';
    end
end

end

