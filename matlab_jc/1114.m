
for n=1:17
    i=n;
    baselineAC=[];
    baselineINA=[];
    for j=1:length(AlldataFinal2(i).exp)
            if AlldataFinal2(i).exp(1).baseline(j)==1 && AlldataFinal2(i).exp(1).acsf(j)==1
                baselineAC=[baselineAC,AlldataFinal2(i).exp(j).selectedpitchcurves(1:1870,:)];
            else if AlldataFinal2(i).exp(1).baseline(j)==1 && AlldataFinal2(i).exp(1).acsf(j)==0
                baselineINA=[baselineINA,AlldataFinal2(i).exp(j).selectedpitchcurves(1:1870,:)];
                end
            end
    end
    
    AlldataFinal2(i).baselineAC=baselineAC;
    AlldataFinal2(i).baselineINA=baselineINA;
end


for i=1:7
    baselineAC=[];
    baselineINA=[];
    for j=1:length(Alldatasecondseven2sigT(i).exp)
    if Alldatasecondseven2sigT(i).exp(1).baseline(j)==1 && Alldatasecondseven2sigT(i).exp(1).acsf(j)==1
        baselineAC=[baselineAC Alldatasecondseven2sigT(i).exp(j).selectedpitchcurves];
    else if Alldatasecondseven2sigT(i).exp(1).baseline(j)==1 && Alldatasecondseven2sigT(i).exp(1).acsf(j)==0
        baselineINA=[baselineINA Alldatasecondseven2sigT(i).exp(j).selectedpitchcurves];
        end
    end
    end
    Alldatasecondseven2sigT(i).baselineAC=baselineAC;
    Alldatasecondseven2sigT(i).baselineINA=baselineINA;
end


figure;plot(std(Alldata1110(4).baselineAC(Alldata1110(4).basevaronset:Alldata1110(4).basevaroffset,:)'),'g')
hold on;plot(std(Alldata1110(4).baselineINA(Alldata1110(4).basevaronset:Alldata1110(4).basevaroffset,:)'),'r')
figure;plot(std(Alldata1110(6).baselineAC(Alldata1110(6).basevaronset:Alldata1110(6).basevaroffset,:)'),'b')
hold on;plot(std(Alldata1110(6).baselineINA(Alldata1110(6).basevaronset:Alldata1110(6).basevaroffset,:)'),'k')
hold on;plot(std(Alldata1110(10).baselineAC(Alldata1110(10).basevaronset:Alldata1110(10).basevaroffset,:)'),'g')
hold on;plot(std(Alldata1110(10).baselineINA(Alldata1110(10).basevaronset:Alldata1110(10).basevaroffset,:)'),'r')
hold on;plot(std(Alldata1110(11).baselineAC(Alldata1110(11).basevaronset:Alldata1110(11).basevaroffset,:)'),'g')
hold on;plot(std(Alldata1110(11).baselineINA(Alldata1110(11).basevaronset:Alldata1110(11).basevaroffset,:)'),'r')
hold on;plot(std(Alldata1110(15).baselineAC(Alldata1110(15).basevaronset:Alldata1110(15).basevaroffset,:)'),'b')
hold on;plot(std(Alldata1110(15).baselineINA(Alldata1110(15).basevaronset:Alldata1110(15).basevaroffset,:)'),'k')
hold on;plot(std(Alldata1110(16).baselineAC(Alldata1110(16).basevaronset:Alldata1110(16).basevaroffset,:)'),'b')
hold on;plot(std(Alldata1110(16).baselineINA(Alldata1110(16).basevaronset:Alldata1110(16).basevaroffset,:)'),'k')

figure;
g=[6 15 16];
for i=1:6
   n=g(i);
    hold on;plot(std(Alldata1110(n).baselineAC(Alldata1110(n).basevaronset:Alldata1110(n).basevaroffset,:)')./std(Alldata1110(n).baselineINA(Alldata1110(n).basevaronset:Alldata1110(n).basevaroffset,:)'),'r')
end


figure;plot(std(Alldata1110(10).baselineAC'))
hold on;plot(std(Alldata1110(10).baselineINA'),'r')

Alldata1110(6).baselineAC=[zeros(1,80) Alldata1110(6).baselineAC];
figure;plot(Alldata1110(4).baselineINA)



for i=1:6
    if i==1 || i==3 || i==4
        lowlimit=2800;
        highlimit=3900;
    else
        lowlimit=2000;
        highlimit=3000;
    end
    Alld2INA(i).exp=jctester2(Alldata_INAtrials(i).exp,lowlimit,highlimit);
end
    



%%%%%% FIGURE 4
% Figure 4A
    slopesACSF=Lesionanaly(Alldata_INAtrials,20,1);
    slopesINA=Lesionanaly(Alldata_INAtrials,20,0);
    figure;plot([1 2],slopesACSF([3 4]),'*')
    hold on;plot([1 2],slopesINA([3 4]),'*','Color','k')
    hold on;plot([3 4 5],slopesACSF([2 5 6]),'*','Color','r')
    hold on;plot([3 4 5],slopesINA([2 5 6]),'*','Color','g')
    ylim([0 0.1])
    xlim([0 6])
    
    
% Figure 4a -- variability reduction is proportional to timescale reduction
% Get timescale data
    for n=1:length(AlldataT)
        i=n;
        amntAC(n)=(mean(std(AlldataT(i).baselineAC(AlldataT(i).basevaronset:AlldataT(i).basevaroffset,:)')));
        amntINA(n)=(mean(std(AlldataT(i).baselineINA(AlldataT(i).basevaronset:AlldataT(i).basevaroffset,:)')));
    end
    figure;plot(amntAC);hold on;plot(amntINA,'r')
    amntratios=(amntINA./amntAC);
    
% Get CV data
    slopesINA=INAanaly(AlldataT,20,0);
    slopesAC=INAanaly(AlldataT,20,1);
    timescaleINA=1./slopesINA;
    timescaleAC=1./slopesAC;
    slopesratios=(timescaleINA./timescaleAC);
    
figure;plot(amntratios([2 4 5 7 8]),slopesratios([2 4 5 7 8]),'*','Color','r') % short
hold on;plot(amntratios([1 3 6 9 10 11]),slopesratios([1 3 6 9 10 11]),'*','Color','k') % long





gg=[];
count=0;
for i=1:313
    if Alldatafirstten2sigT(4).baselineAC(350,i)<3700
        count=count+1;
        gg(:,count)=Alldatafirstten2sigT(4).baselineAC(:,i);
    end
end


%xcorrelations

for i=1:12
    gg=Alldata2(1).ind_longnotes(i);
    pitchAC=Alldata2(Alldata2(gg).baselineAC(Alldata2(gg).basevaronset+40:Alldata2(gg).basevaroffset-40,:));
    pitchINA=Alldata2(Alldata2(gg).baselineINA(Alldata2(gg).basevaronset+40:Alldata2(gg).basevaroffset-40,:));
    figure;plot(xcorr(pitchAC))
    hold on;plot(xcorr(pitchINA),'r')
end
