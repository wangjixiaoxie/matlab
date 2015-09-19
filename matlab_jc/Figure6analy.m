Figure6analy
%

g=[1];
for n=1:length(g)
    i=g(n);
    AD2final(i).exp=jctester2(Alldatasecondseven1sig(i).exp,2900,3900);
end

Alldatafirstten2sigT=Alldatafirstten1sig;
for i=1:10
    for j=1:length(Alldatafirstten2sigT(i).exp)
        Alldatafirstten2sigT(i).exp(j).rawpitchcurves=Alldatafirstten2sig(i).exp(j).pitchcurves;
    end
end
for k=1:10
Alldatafirstten2sigT(k).exp=jctester3(Alldatafirstten2sigT(k).exp,Alldatafirstten2sigT(k).basevaronset, Alldatafirstten2sigT(k).basevaroffset);
end