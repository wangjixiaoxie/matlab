function Figure4analy(Alldata,targregion,conting)

%%%%% This code tests to see which experiments have good baseline
%%%%% inactivations
% n=12;
% alls=find(Alldata(n).exp(1).baseline==0);
% figure;hold on;
% for i=1:alls(1)-1 % for all baseline groups
%     if Alldata(n).exp(1).acsf(i)==0
%         % Plot CVtrace
%         plot(std(Alldata(n).exp(i).selectedpitchcurves')./mean(Alldata(n).exp(i).selectedpitchcurves'),'r')
%     else
%         plot(std(Alldata(n).exp(i).selectedpitchcurves')./mean(Alldata(n).exp(i).selectedpitchcurves'),'k')
%     end
% end
ACSF=0;
conting=20;
targregion=3;
Alld=Alldata([4 6 10 11 15 16]);

[longer,middle,allcontinga]=Fig3ContingSimulator(Alld,targregion,conting,ACSF);