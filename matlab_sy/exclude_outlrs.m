function [fvnew]=exclude_outlrs(fv, bds);

    vls=getvals(fv,1,'trig')
    stdv=std(vls(:,2));
    mnvl=mean(vls(:,2));
    diff=vls(:,2)-mnvl;
    ind=find(abs(diff)<2*stdv);
    fvnew=fv(ind);
    


