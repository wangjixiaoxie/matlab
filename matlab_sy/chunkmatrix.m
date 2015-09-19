
function [chkarray]=chunkmatrix(vals,stimes,etimes,tcl)
    
    
    for ii=1:length(stimes)
        ind=find(vals(:,tcl)>=stimes(ii)&vals(:,tcl)<=etimes(ii))
        chkarray{ii}=vals(ind,:)
    end
        


