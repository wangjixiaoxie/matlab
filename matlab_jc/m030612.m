TW data - already processed
load /bulbul1/r87g80rvsd.mat
for i=1:length(vls)
vls(i,6)=fvpt_com(i).TRIG(1);
end

CTind=crct_com(1:53246);
FBind=crfb_com(1:8919);
vlsCT=vls(CTind,:);
vlsFB=vls(FBind,:);
vlsCTexp1=vlsCT(find(vlsCT(:,1)<734500),:);
vlsFBexp1=vlsFB(find(vlsFB(:,1)<734500),:);
vlsCTexp2=vlsCT(find(vlsCT(:,1)>734500),:);
vlsFBexp2=vlsFB(find(vlsFB(:,1)>734500),:);
[IntExp(5).deltadayFBlearning,IntExp(5).deltanightFBlearning,IntExp(5).AllDaysFBlearning,IntExp(5).FFmorningfitFBlearning,IntExp(5).FFeveningfitFBlearning]...
    =daynightTW(vlsFBexp1);
[IntExp(6).deltadayFBlearning,IntExp(6).deltanightFBlearning,IntExp(6).AllDaysFBlearning,IntExp(6).FFmorningfitFBlearning,IntExp(6).FFeveningfitFBlearning]...
    =daynightTW(vlsFBexp2);
[IntExp(5).deltadayCTlearning,IntExp(5).deltanightCTlearning,IntExp(5).AllDaysCTlearning,IntExp(5).FFmorningfitCTlearning,IntExp(5).FFeveningfitCTlearning]...
    =daynightTW(vlsCTexp1);
[IntExp(6).deltadayCTlearning,IntExp(6).deltanightCTlearning,IntExp(6).AllDaysCTlearning,IntExp(6).FFmorningfitCTlearning,IntExp(6).FFeveningfitCTlearning]...
    =daynightTW(vlsCTexp2);

IntExp(7).deltadayFBlearning=IntExp(5).deltadayFBlearning(13:14);IntExp(7).deltadayFBbaseline=IntExp(5).deltadayFBlearning(11:2);
IntExp(7).deltanightFBlearning=IntExp(5).deltanightFBlearning(13:14);IntExp(7).deltanightFBbaseline=IntExp(5).deltanightFBlearning(11:12);
IntExp(7).FFmorningfitFBlearning=IntExp(5).FFmorningfitFBlearning(13:14);IntExp(7).FFmorningfitFBbaseline=IntExp(5).FFmorningfitFBlearning(11:12);
IntExp(7).AllDaysFBlearning=IntExp(5).AllDaysFBlearning(13:14);IntExp(7).AllDaysFBbaseline=IntExp(5).AllDaysFBlearning(11:12);
IntExp(7).FFeveningfitFBlearning=IntExp(5).FFeveningfitFBlearning(13:14);IntExp(7).FFeveningfitFBbaseline=IntExp(5).FFeveningfitFBlearning(11:12);
IntExp(7).deltadayCTlearning=IntExp(5).deltadayCTlearning(13:14);IntExp(7).deltadayCTbaseline=IntExp(5).deltadayCTlearning(11:12);
IntExp(7).deltanightCTlearning=IntExp(5).deltanightCTlearning(13:14);IntExp(7).deltanightCTbaseline=IntExp(5).deltanightCTlearning(11:12);
IntExp(7).FFmorningfitCTlearning=IntExp(5).FFmorningfitCTlearning(13:14);IntExp(7).FFmorningfitCTbaseline=IntExp(5).FFmorningfitCTlearning(11:12);
IntExp(7).AllDaysCTlearning=IntExp(5).AllDaysCTlearning(13:14);IntExp(7).AllDaysCTbaseline=IntExp(5).AllDaysCTlearning(11:12);
IntExp(7).FFeveningfitCTlearning=IntExp(5).FFeveningfitCTlearning(13:14);IntExp(7).FFeveningfitCTbaseline=IntExp(5).FFeveningfitCTlearning(11:12);

IntExp(5).deltadayFBlearning=IntExp(5).deltadayFBlearning(4:10);IntExp(5).deltadayFBbaseline=IntExp(5).deltadayFBlearning(1:2);
IntExp(5).deltanightFBlearning=IntExp(5).deltanightFBlearning(4:10);IntExp(5).deltanightFBbaseline=IntExp(5).deltanightFBlearning(1:2);
IntExp(5).FFmorningfitFBlearning=IntExp(5).FFmorningfitFBlearning(4:10);IntExp(5).FFmorningfitFBbaseline=IntExp(5).FFmorningfitFBlearning(1:2);
IntExp(5).AllDaysFBlearning=IntExp(5).AllDaysFBlearning(4:10);IntExp(5).AllDaysFBbaseline=IntExp(5).AllDaysFBlearning(1:2);
IntExp(5).FFeveningfitFBlearning=IntExp(5).FFeveningfitFBlearning(4:10);IntExp(5).FFeveningfitFBbaseline=IntExp(5).FFeveningfitFBlearning(1:2);
IntExp(5).deltadayCTlearning=IntExp(5).deltadayCTlearning(4:10);IntExp(5).deltadayCTbaseline=IntExp(5).deltadayCTlearning(1:2);
IntExp(5).deltanightCTlearning=IntExp(5).deltanightCTlearning(4:10);IntExp(5).deltanightCTbaseline=IntExp(5).deltanightCTlearning(1:2);
IntExp(5).FFmorningfitCTlearning=IntExp(5).FFmorningfitCTlearning(4:10);IntExp(5).FFmorningfitCTbaseline=IntExp(5).FFmorningfitCTlearning(1:2);
IntExp(5).AllDaysCTlearning=IntExp(5).AllDaysCTlearning(4:10);IntExp(5).AllDaysCTbaseline=IntExp(5).AllDaysCTlearning(1:2);
IntExp(5).FFeveningfitCTlearning=IntExp(5).FFeveningfitCTlearning(4:10);IntExp(5).FFeveningfitCTbaseline=IntExp(5).FFeveningfitCTlearning(1:2);

IntExp(8).deltadayFBlearning=IntExp(6).deltadayFBlearning(16:19);IntExp(8).deltadayFBbaseline=IntExp(6).deltadayFBlearning(14:15);
IntExp(8).deltanightFBlearning=IntExp(6).deltanightFBlearning(16:19);IntExp(8).deltanightFBbaseline=IntExp(6).deltanightFBlearning(14:15);
IntExp(8).FFmorningfitFBlearning=IntExp(6).FFmorningfitFBlearning(16:19);IntExp(8).FFmorningfitFBbaseline=IntExp(6).FFmorningfitFBlearning(14:15);
IntExp(8).AllDaysFBlearning=IntExp(6).AllDaysFBlearning(16:19);IntExp(8).AllDaysFBbaseline=IntExp(6).AllDaysFBlearning(14:15);
IntExp(8).FFeveningfitFBlearning=IntExp(6).FFeveningfitFBlearning(16:19);IntExp(8).FFeveningfitFBbaseline=IntExp(6).FFeveningfitFBlearning(14:15);
IntExp(8).deltadayCTlearning=IntExp(6).deltadayCTlearning(16:19);IntExp(8).deltadayCTbaseline=IntExp(6).deltadayCTlearning(14:15);
IntExp(8).deltanightCTlearning=IntExp(6).deltanightCTlearning(16:19);IntExp(8).deltanightCTbaseline=IntExp(6).deltanightCTlearning(14:15);
IntExp(8).FFmorningfitCTlearning=IntExp(6).FFmorningfitCTlearning(16:19);IntExp(8).FFmorningfitCTbaseline=IntExp(6).FFmorningfitCTlearning(14:15);
IntExp(8).AllDaysCTlearning=IntExp(6).AllDaysCTlearning(16:19);IntExp(8).AllDaysCTbaseline=IntExp(6).AllDaysCTlearning(14:15);
IntExp(8).FFeveningfitCTlearning=IntExp(6).FFeveningfitCTlearning(16:19);IntExp(8).FFeveningfitCTbaseline=IntExp(6).FFeveningfitCTlearning(14:15);


IntExp(6).deltadayFBlearning=IntExp(6).deltadayFBlearning(3:end);IntExp(6).deltadayFBbaseline=IntExp(6).deltadayFBlearning(1:2);
IntExp(6).deltanightFBlearning=IntExp(6).deltanightFBlearning(3:end);IntExp(6).deltanightFBbaseline=IntExp(6).deltanightFBlearning(1:2);
IntExp(6).FFmorningfitFBlearning=IntExp(6).FFmorningfitFBlearning(3:end);IntExp(6).FFmorningfitFBbaseline=IntExp(6).FFmorningfitFBlearning(1:2);
IntExp(6).AllDaysFBlearning=IntExp(6).AllDaysFBlearning(3:end);IntExp(6).AllDaysFBbaseline=IntExp(6).AllDaysFBlearning(1:2);
IntExp(6).FFeveningfitFBlearning=IntExp(6).FFeveningfitFBlearning(3:end);IntExp(6).FFeveningfitFBbaseline=IntExp(6).FFeveningfitFBlearning(1:2);
IntExp(6).deltadayCTlearning=IntExp(6).deltadayCTlearning(3:end);IntExp(6).deltadayCTbaseline=IntExp(6).deltadayCTlearning(1:2);
IntExp(6).deltanightCTlearning=IntExp(6).deltanightCTlearning(3:end);IntExp(6).deltanightCTbaseline=IntExp(6).deltanightCTlearning(1:2);
IntExp(6).FFmorningfitCTlearning=IntExp(6).FFmorningfitCTlearning(3:end);IntExp(6).FFmorningfitCTbaseline=IntExp(6).FFmorningfitCTlearning(1:2);
IntExp(6).AllDaysCTlearning=IntExp(6).AllDaysCTlearning(3:end);IntExp(6).AllDaysCTbaseline=IntExp(6).AllDaysCTlearning(1:2);
IntExp(6).FFeveningfitCTlearning=IntExp(6).FFeveningfitCTlearning(3:end);IntExp(6).FFeveningfitCTbaseline=IntExp(6).FFeveningfitCTlearning(1:2);

