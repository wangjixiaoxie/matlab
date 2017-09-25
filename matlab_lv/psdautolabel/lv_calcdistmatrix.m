function [distmtx, normby] = lv_calcdistmatrix(basismtx,allpsds,normby)

basissylnum=size(basismtx,1);
distmtx=zeros(size(allpsds,1),basissylnum);
% distmtx = allpsds * basismtx';

for i=1:(size(allpsds,1))
    for j=1:(basissylnum)
        distmtx(i,j) = max(xcorr(allpsds(i,:),basismtx(j,:),15));
        
        %%sum(sqrt(((allpsds(i,:)-basismtx(j,:)).^ 2)));
    end
end

% distmtx = 1-(distmtx./max(max(distmtx))); %convert to similarity matrix  
if normby == 0
    normby = max(max(distmtx));
end
distmtx = distmtx./normby;
