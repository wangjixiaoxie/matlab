function [medcon]=jc_consensus(f_est,q)
%Run this program after running ifdv with a break point around line 99

a=size(f_est);

% Remove all points that aren't in lock from f_est
for i=1:a(1)
    for j=1:a(2)
        if q(i,j)==0
            f_est(i,j)=0;
        end
    end
end

consensus(1)=1; %initialize consensus

%For all adjacent frequency bins that are in lock, determine the difference
%between the estimated frequency in each bin.
for j=1:a(2) % for each time bin (i.e. each column)
    for i=1:a(1)-1
        if f_est(i,j)~=0    %If this point is locked and...
            if f_est(i+1,j)~=0  %If the next point is locked, then...
                cons=1/abs(f_est(i,j)-f_est(i+1,j));
                consensus=[consensus cons];
            end
        end
    end
end
consensus=consensus(2:length(consensus)); %remove the initialization point
medcon=median(consensus); %Take the median value
                
            
                