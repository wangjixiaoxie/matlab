%output sortvals is a cell array, this gets called by chunkdata 2

function [sortvals]=sort_times(vals,timevec)
for ii=1:length(timevec)
    ind=find(vals(:,1)>timevec(ii,1)&vals(:,1)<timevec(ii,2))
    sortvals{ii}=vals(ind,:);
end