DATA SET SUMMARY
%   BF Inactivations that decrease CV- Long (6 15 16); Short (10 11)
        % also bk61w42
%   BF Notes in which FD decreases CV - Long (3); Short (2,4)
%   ZF Notes in which (Lesions, FD) decrease CV - (1 through 6) (7 is my bird)

%%%%% BF inactivations
    Alldata2=[Alldatafirstten2sigT Alldatasecondseven2sigT];
        % The same dataset with pitchcurves calculated at 1sigma is also there
    % Alldata2 contains BF inactivation data
    % toffsets contains targeting distributions for Alldata2
    
    % Goodinactivations contains successful (deltaCV decreases) baseline inactivation runs
        % Goodinactivations.Allinact selected as birds any decrease in CV
        % The following 5 were selected as trials that show large CV decreases.
            % 1-bk61w42 (Allinact(1) --- AC10vs.INA5)
            % 2-bk61w42 (Allinact(1) --- AC5vs.INA3)
            % 3-Alldata2(6) (Allinact(2) --- AC7vs.INA5)
            % 4-Alldata2(16) (Allinact(4) --- AC3vs.INA2)
            % 5-Alldata2(11) (Allinact(6) --- AC1vs.INA1)
%%%%% BF lesions  
    % AlldataBFlesion
%%%%% ZF lesions (and #7 is non-lesioned - my bird)
    % AlldataZFlesion
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DATA ANALYSIS
                % 1209.m and 1114.m are good for answering any questions
                % left unanswering by the code below -- 1209.m is highly
                % specific

    g=[Alldata2(1).ind_allnotes 17]; % notes that you're testing

% Actual learning (i.e. shifted pitch curves) in three raw-->normalized formats
    n=4; % the learning block experiment you're testing
    [actshift,normshift,pcnormshift,maxshiftlocus]=jc_actualshift(Alldata2(n),1,0.8,24);
% Predicted learning
    avf=ContingSim(toffsets(6).data,Alldata2(6).baselineAC,70);
% Predicting position of maximal shift from mean of targeting distribution
    mediantargeting=targetingdistn(Alldata2(g),1,0.8,24,0); 
    for i=1:length(g)
        [a,b,c,maxshiftlocus(i)]=jc_actualshift(Alldata2(i),1,0.8,24);
    end
    figure;plot(mediantargeting,maxshiftlocus,'*')
    
% Autocorrelation based analysis/quantification
            % Edit g in these programs to only get the desired experiments
        % These programs plot the datapoints that they return.
            % Edit programs internally to plot to autocorrelation data.
    [slopesAC,slopesINA,peaksAC,peaksINA]=xcorranalyINA(Alldata2);
    [slopesAC,slopesINA,peaksAC,peaksINA]=xcorranalyLESION(AlldataZFlesion(1:6));
    
% Simulate learning at three distinct points (beginning, middle, end) using ALL variation present
    figure;LMANanalyPlot(Alldata2,20,1); % Change program internally to do specific notes - e.g. long vs. short
        % 1 means ACSF, 0 means INA, 2 means UDpre, 3 means Dpre, 4 means UDpost
        
    % Compare slopes at the three different points.
        [slopebegin,slopemiddle,slopeend]=LMANanalyPlot2(Alldata2(g),20,1);
        
% Simulate learning at targeting distn using ALL variation present
    [a,targetdist]=targetingdistn(Alldata2(g),1,0.8,24,0);

% Decrease in CV vs. Decrease in Timescale PLOT
    CVTimescaleplot(Alldata2); % - the 4th figure is the one we want
    
% Plotting the spectrograms and overlaying the shifts and targeting distn
    Figure5bk50.m --- has commands but not intended to run
    Figure5pk37.m --- has commands but not intended to run
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compare different methods for obtaining pitch curves
    % 1 sigma, 2 sigma --- see above
    jc_hilbert(shifted,2000,3000,2); 
        % takes the hilbert transform butter bandpass filtering 2000 to
        % 3000 Hz and then butter lowpass smoothing the end result to
        % eliminate fluctuations faster than 2ms.
% Compare them using autocorrelation or simulation method

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%








% Simulate learning at targeting distn using only LMAN variation
--- CODE BELOW ATTEMPTS TO DO THIS --- need to use good targeting distn and make other modifications
    1215.m, 1215b.m,1212.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% Problem: how to account for offsets when they occur before or after the point
% in the note where you start calculating the LMAN contribution --- perhaps
% just have the distribution (standev) w/o the location.


% LMAN contribution pitch curves
LMAN(1).data=jc_LMANvariation(ACSF(1).residSTFT(450:730,:),INA(1).residSTFT(450:730,:))'; % bk61w42 
LMAN(2).data=jc_LMANvariation(ACSF(2).residSTFT(450:730,:),INA(2).residSTFT(450:730,:))'; % bk61w42
LMAN(3).data=jc_LMANvariation(ACSF(3).residSTFT(270:580,:),INA(3).residSTFT(270:580,:))'; % Alldata2(6)
LMAN(4).data=jc_LMANvariation(ACSF(4).residSTFT(730:1070,:),INA(4).residSTFT(730:1070,:))'; % Alldata2(16)
LMAN(5).data=jc_LMANvariation(ACSF(5).residSTFT(280:400,:),INA(5).residSTFT(280:400,:))'; % Alldata2(11)

% All variation pitch curves
ALL(1).data=ACSF(1).pitchSTFT(450:730,:);
ALL(2).data=ACSF(2).pitchSTFT(450:730,:);
ALL(3).data=Alldata2(6).baselineAC(270:580,:);
ALL(4).data=Alldata2(16).baselineAC(730:1070,:);
ALL(5).data=Alldata2(11).baselineAC(280:400,:);

Simtoffs(1).data=ones(1,50)*100;
Simtoffs(2).data=Simtoffs(1).data;
Simtoffs(3).data=Simtoffs(1).data;
Simtoffs(4).data=Simtoffs(1).data;
Simtoffs(5).data=Simtoffs(1).data;
% Simtoffs(3).data=toffsets(6).data-270;
% Simtoffs(3).data(find(Simtoffs(3).data<33))=33;
% Simtoffs(3).data(find(Simtoffs(3).data>size(ALL(3).data,1)))=size(ALL(3).data,1);
% Simtoffs(4).data=toffsets(16).data-730;
% Simtoffs(4).data(find(Simtoffs(4).data<33))=33;
% Simtoffs(4).data(find(Simtoffs(4).data>size(ALL(4).data,1)))=size(ALL(4).data,1);
% Simtoffs(5).data=toffsets(11).data-280;
% Simtoffs(5).data(find(Simtoffs(5).data<33))=33;
% Simtoffs(5).data(find(Simtoffs(5).data>size(ALL(5).data,1)))=size(ALL(5).data,1);

Simtoffs(3).data=toffsets(6).data;
Simtoffs(3).data(find(Simtoffs(3).data<33))=33;
Simtoffs(3).data(find(Simtoffs(3).data>size(ALL(3).data,1)))=size(ALL(3).data,1);
Simtoffs(4).data=toffsets(16).data;
Simtoffs(4).data(find(Simtoffs(4).data<33))=33;
Simtoffs(4).data(find(Simtoffs(4).data>size(ALL(4).data,1)))=size(ALL(4).data,1);
Simtoffs(5).data=toffsets(11).data;
Simtoffs(5).data(find(Simtoffs(5).data<33))=33;
Simtoffs(5).data(find(Simtoffs(5).data>size(ALL(5).data,1)))=size(ALL(5).data,1);



% Simulate both

for i=1:5
    avALL(i).data=ContingSim2(Simtoffs(i).data,ALL(i).data,70);
    avLMAN(i).data=ContingSim2(Simtoffs(i).data,LMAN(i).data,70);
    figure;plot(avALL(i).data./max(avALL(i).data))
    hold on;plot(avLMAN(i).data./max(avLMAN(i).data),'r')
end
% Convert from curves into one-sided doohickeys
% Plot *actual* shifts
targetingplotfigure2