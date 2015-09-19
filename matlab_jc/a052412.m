load /bulbul1/IntExp030612.mat
for i=1:4
    sdCT(i)=std(IntExp(i).FFdataCTbaseline);
    mnCT(i)=mean(IntExp(i).FFdataCTbaseline);
    sdFB(i)=std(IntExp(i).FFdataFBbaseline);
    mnFB(i)=mean(IntExp(i).FFdataFBbaseline);
end
for i=5:8
    FvalsCT=[];
    FvalsFB=[];
    for j=1:length(IntExp(i).AllDaysFBbaseline)
        FvalsCT=[FvalsCT;IntExp(i).AllDaysCTbaseline(j).FFdata];
        FvalsFB=[FvalsFB;IntExp(i).AllDaysFBbaseline(j).FFdata];
    end
    sdCT(i)=std(FvalsCT);
    mnCT(i)=mean(FvalsCT);
    sdFB(i)=std(FvalsFB);
    mnFB(i)=mean(FvalsFB);
end
sdFB=[sdFB(1:4) mean(sdFB(5:8))];
sdCT=[sdCT(1:4) mean(sdCT(5:8))];

mean(sdFB./sdCT) % 0.95
mean(sdFB./mnFB)./mean(sdCT./mnCT) % 0.9243

