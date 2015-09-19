function Alldata=baselinecat(Alldata)
for i=1:length(Alldata)
    CurrExp=Alldata(i).exp;
    baselineAC=[];
    baselineINA=[];
    n=CurrExp(1).baseline(1);
    j=1;
    sizer=size(CurrExp(1).selectedpitchcurves,1);
    while n==1
        if CurrExp(1).acsf(j)==1
            baselineAC=[baselineAC CurrExp(j).selectedpitchcurves(1:sizer,:)];
        else
            baselineINA=[baselineINA CurrExp(j).selectedpitchcurves(1:sizer,:)];
        end
        j=j+1;
        n=CurrExp(1).baseline(j);
    end
    Alldata(i).baselineAC=baselineAC;
    Alldata(i).baselineINA=baselineINA;
end
