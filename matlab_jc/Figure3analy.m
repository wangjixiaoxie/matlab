function [slopes,aligneddata]=Figure3analy(targregion,Alldata,cont,ACSF)

% Calls Fig3ContingSimulator.m
% Calls Fig3Aligner.m



% Parameters
    % targregion: 1=beginning, 2=middle, 3=end
    % cont: e.g. 20 - percentile from mean to exceed 20=mean (30,70)
    % Alldata - selected experiments e.g. Alldata(Alldata(1).ind_longnotes)
    % ACSF - 1 if ACSF, 0 if inactivations


% Select the curves that meet the above and below contingencies
for mm=1:2
    if mm==1
        conting=50+cont;
        [longer,middle,allcontinga]=Fig3ContingSimulator(Alldata,targregion,conting,ACSF,mm);
    else
        conting=50-cont;
        [longer,middle,allcontingb]=Fig3ContingSimulator(Alldata,targregion,conting,ACSF,mm);
    end
end

% Take the mean of above and below
allcontingcurves=(allcontinga+allcontingb)/2;

% Align the notes
[tracee]=Fig3Aligner(Alldata,allcontingcurves,targregion,longer,middle);

% Calculate the slopes of decay
x=(0:1/8:119/8);
for i=1:length(Alldata)
    kk1=polyfit(x,tracee(i,1:120),1);
    slopes(i)=kk1(1);
end



% Plot mean and standard error above and below
for i=1:1000
    xax(i)=i/8;
end
% %%% Only for ACSF/non-ACSF runs
% if ACSF==1
%      plot(xax,mean(tracee),'k')
% else
%      plot(xax,mean(tracee),'r')
% end
aligneddata=tracee;


% Figure 3 plotting commands
llong=120;
lshort=65;
     figure;plot(xax(1:llong),mean(tracee(Alldata(1).ind_longnotes,1:llong)),'b');
    hold on;plot(xax(1:llong),mean(tracee(Alldata(1).ind_longnotes,1:llong))+std(tracee(Alldata(1).ind_longnotes,1:llong))/sqrt(length(Alldata(1).ind_longnotes)),'b')
    hold on;plot(xax(1:llong),mean(tracee(Alldata(1).ind_longnotes,1:llong))-std(tracee(Alldata(1).ind_longnotes,1:llong))/sqrt(length(Alldata(1).ind_longnotes)),'b')
    hold on;plot(xax(1:lshort),mean(tracee(Alldata(1).ind_shortnotes,1:lshort)),'r');
    hold on;plot(xax(1:lshort),mean(tracee(Alldata(1).ind_shortnotes,1:lshort))+std(tracee(Alldata(1).ind_shortnotes,1:lshort))/sqrt(length(Alldata(1).ind_shortnotes)),'r')
    hold on;plot(xax(1:lshort),mean(tracee(Alldata(1).ind_shortnotes,1:lshort))-std(tracee(Alldata(1).ind_shortnotes,1:lshort))/sqrt(length(Alldata(1).ind_shortnotes)),'r')
ylim([0 1.05])

% Figure 4 plotting commands
llong=120;
lshort=65;
indlong=[1 4 5];
indshort=[2 3];
if ACSF==0
    cL='k';
    cS='r';
else
    cL='b';
    cS='g';
end
    hold on;plot(xax(1:llong),mean(tracee(indlong,1:llong)),cL)
    hold on;plot(xax(1:llong),mean(tracee(indlong,1:llong))+std(tracee(indlong,1:llong))/sqrt(length(Alldata)),cL)
    hold on;plot(xax(1:llong),mean(tracee(indlong,1:llong))-std(tracee(indlong,1:llong))/sqrt(length(Alldata)),cL)
    hold on;plot(xax(1:lshort),mean(tracee(indshort,1:lshort)),cS);
    hold on;plot(xax(1:lshort),mean(tracee(indshort,1:lshort))+std(tracee(indshort,1:lshort))/sqrt(length(Alldata)),cS)
    hold on;plot(xax(1:lshort),mean(tracee(indshort,1:lshort))-std(tracee(indshort,1:lshort))/sqrt(length(Alldata)),cS)
    ylim([0 1.05])

