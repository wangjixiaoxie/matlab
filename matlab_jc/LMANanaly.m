function [slopes]=LMANanaly(Alldata,cont,ACSF)

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
x=(0:1/8:59/8);
for i=1:size(traces,1)
    kk1=polyfit(x,traces(i,1:60),1);
    slopes(i)=-1*kk1(1);
end
