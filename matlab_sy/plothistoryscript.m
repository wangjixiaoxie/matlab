
figure
mnoutinds=[1 2 4 5]
colinds={'k' 'r' 'c' 'b'}

for ii=1:length(mnoutinds)
    crmnout=mnout(mnoutinds(ii));
    for inds=1:length(crmnout.acsh)
        shdiff=abs(crmnout.acsh{inds})-abs(crmnout.mush{inds});
        rvdiff=abs(crmnout.murev{inds})-abs(crmnout.acrev{inds});
    
    plot((shdiff-rvdiff)/2,(shdiff+rvdiff)/2,'ko','MarkerSize',5,'MarkerFaceColor','k');
    hold on;
    end
    end
    
