function [alignedtrace]=Fig3Aligner(Alldata,allcontingcurves,targregion,longer,middle)

if targregion==1
    tracee=zeros(length(Alldata),1000);
    for aa=1:length(Alldata);
        for i=1:longer(aa)   
            tracee(aa,i)=(allcontingcurves(aa,i));
        end
    end
end
if targregion==2
    tracee=zeros(length(Alldata),1000);
    for aa=1:length(Alldata);
        for i=1:middle(aa)-1  
            tracee(aa,i)=(allcontingcurves(aa,middle(aa)-i)+allcontingcurves(aa,middle(aa)+i))/2;
        end
    end
end
if targregion==3
    tracee=zeros(length(Alldata),1000);
    for aa=1:length(Alldata);
        for i=1:longer(aa)-1
            tracee(aa,i)=(allcontingcurves(aa,longer(aa)-(i)));
        end
    end
end
alignedtrace=tracee;