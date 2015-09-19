%9/15/08 purpose of this function is to sort all the muscimol runs into
%following categories-
%baseline runs, as defined by being 1 standard deviation from the mean
%maxruns, within 1.5 standard deviation of maximum, combined 1.5 standard
%deviation run.
function [avls]=sortmuruns(avls);
   asympfac=1; 
%loop through notes
for nt_nm=1:length(avls.muanal)
   
    initmean=avls.initmean{nt_nm};
    initsd=avls.initsd{nt_nm};
    numepochs=length(avls.muanal{nt_nm})
    for ep_nm=1:numepochs
        
        muind=avls.muanal{nt_nm}{ep_nm}
        muvls=avls.acmean(nt_nm,muind);
        
        %find the values greater than onesd.
        diff=abs(muvls-initmean);
        basind=find(diff<initsd);
        shiftind=find(diff>=initsd);
        avls.basind{nt_nm}{ep_nm}=muind(basind);
        avls.shiftind{nt_nm}{ep_nm}=muind(shiftind);
        %now find the values around the maximum shift
        max_shift=max(diff);
        diff_mx=abs(max_shift-diff);
        asympind=find(diff_mx<=asympfac*initsd&diff>=initsd);
        avls.asympind{nt_nm}{ep_nm}=muind(asympind);
    end
end
end

        %check whether mumean is within one standard deviation of the
        %original mean
        
        
        
        %check whether mumean is within 1.5 standard deviations of the
        %maximum mumean.