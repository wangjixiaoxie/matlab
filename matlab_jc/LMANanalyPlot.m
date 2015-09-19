function [slopes]=LMANanalyPlot(Alldata,cont,ACSF)
% Simulation of learning at three points using all variation and ignoring
% targeting imprecision

% Calls Fig3ContingSimulator.m
% Calls Fig3Aligner.m

% Parameters
    % targregion: 1=beginning, 2=middle, 3=end
    % cont: e.g. 20 - percentile from mean to exceed 20=mean (30,70)
    % Alldata - selected experiments e.g. Alldata(Alldata(1).ind_longnotes)
    
    % ACSF - 1 if ACSF, 0 if inactivations
    % ACSF - 2 if UDpre, 3 if Dpre, 4 if UDpost --- code in Fig3ContingSimulator.m
    
% For beginning, middle, and end
for targregion=1:3
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
    [tracee(targregion).traces]=Fig3Aligner(Alldata,allcontingcurves,targregion,longer,middle);
end

% Weight equally each place in the note where I'm measuring decay
for i=1:size(tracee(1).traces,1)
    traces(i,:)=(tracee(1).traces(i,:)+2*tracee(2).traces(i,:)+tracee(3).traces(i,:))./4;
end


% Calculate the slopes of decay - first 5ms 
x=(0:1/8:39/8);
for i=1:size(traces,1)
    kk1=polyfit(x,traces(i,1:40),1);
    slopes(i)=-1*kk1(1);
end

for i=1:1000
    xax(i)=i/8;
end

%%%% PLOT
tracee=traces;
% Figure 4 plotting commands
llong=100;
lshort=65;
indlong=[1 6 7 8 9 12 13 14 15 16 17]; % 1 4 5 for INAs, 1 3 for BFLs
indshort=[4 5 10 11];
%indshort=Alldata(1).ind_shortnotes; % 2 3 for INAs, 2 4 for BFLs
if ACSF==1 %1
    cL='k';
    cS='g';
else
    cL='b';
    cS='r';
end
    hold on;plot(xax(1:llong),mean(tracee(indlong,1:llong)),cL)
    hold on;plot(xax(1:llong),mean(tracee(indlong,1:llong))+std(tracee(indlong,1:llong))/sqrt(length(Alldata)),cL)
    hold on;plot(xax(1:llong),mean(tracee(indlong,1:llong))-std(tracee(indlong,1:llong))/sqrt(length(Alldata)),cL)
    hold on;plot(xax(1:lshort),mean(tracee(indshort,1:lshort)),cS);
    hold on;plot(xax(1:lshort),mean(tracee(indshort,1:lshort))+std(tracee(indshort,1:lshort))/sqrt(length(Alldata)),cS)
    hold on;plot(xax(1:lshort),mean(tracee(indshort,1:lshort))-std(tracee(indshort,1:lshort))/sqrt(length(Alldata)),cS)
    ylim([0 1.05])
