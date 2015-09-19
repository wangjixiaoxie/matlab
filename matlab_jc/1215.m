
n=6;
exp=11;
experiment=Alldata2(exp).exp;
countina=0; countac=0;
for i=1:length(experiment)
    if experiment(1).baseline(i)==1 && experiment(1).acsf(i)==1
        countac=countac+1;
        Allinac(n).ACSF(countac).rawdata=experiment(i).rawdata;
        Allinac(n).ACSF(countac).rawpitchcurves=experiment(i).rawpitchcurves;
        Allinac(n).ACSF(countac).selectedpitchcurves=experiment(i).selectedpitchcurves;
    else
        if experiment(1).baseline(i)==1 && experiment(1).acsf(i)==0
            countina=countina+1;
            Allinac(n).INA(countina).rawdata=experiment(i).rawdata;
            Allinac(n).INA(countina).rawpitchcurves=experiment(i).rawpitchcurves;
            Allinac(n).INA(countina).selectedpitchcurves=experiment(i).selectedpitchcurves;
        end
    end
end



for i=1
for j=5:length(Allinact(i).ACSF)
Allinact(i).ACSF(j).residuals=jc_residuals(Allinact(i).ACSF(j).selectedpitchcurves);
end
for k=3:length(Allinact(i).INA)
Allinact(i).INA(k).residuals=jc_residuals(Allinact(i).INA(k).selectedpitchcurves);
end
end