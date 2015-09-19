%% MANUAL - TO compile multiple batch files into one structure, do following:
% 1) first run heat map function for all data. After run, take output
% strcutures and compile by running script below:

% INPUTS - what to call this data.
% fieldname='delay0';
% load DATSTRUCT;
% load ../Params.mat;

FracTrialsForCatch=0.5;

BirdDir='/bluejay3/lucas/birds/wh73pk61/';

DirsOfInterest={'030115_HVCChR2_HVCStim_StimON_quick/lt_Opto_Stim_analy_batch.labeled.all.delay140_10Mar2015_2302/PLOT_StimCatch_StimNotCatch_10Mar2015_2305/TimeWindow_10Mar2015_2308/HEATMAP_11Mar2015_1802',...
    '030115_HVCChR2_HVCStim_StimON_quick/lt_Opto_Stim_analy_batch.labeled.all.delay200_02Mar2015_1059/PLOT_StimCatch_StimNotCatch_02Mar2015_1101/TimeWindow_10Mar2015_1307/HEATMAP_11Mar2015_1804',...
    '030115_HVCChR2_HVCStim_StimON_quick/lt_Opto_Stim_analy_batch.labeled.all.delay220_02Mar2015_1110/PLOT_StimCatch_StimNotCatch_02Mar2015_1114/TimeWindow_09Mar2015_1423/HEATMAP_11Mar2015_1806',...
    '030115_HVCChR2_HVCStim_StimON_quick/lt_Opto_Stim_analy_batch.labeled.all.delay80_02Mar2015_1154/PLOT_StimCatch_StimNotCatch_02Mar2015_1158/TimeWindow_10Mar2015_1244/HEATMAP_11Mar2015_1936',...
    '030215_HVCChR2_HVCStim_StimON_quick/lt_Opto_Stim_analy_batch.labeled.all.nodelay_02Mar2015_1214/PLOT_StimCatch_StimNotCatch_02Mar2015_1216/TimeWindow_10Mar2015_1253/HEATMAP_11Mar2015_1938'};

    
    
% 
DATSTRUCTComp=struct;
PARAMSComp=struct;

%% LOAD STUFF

NumStructs=length(DirsOfInterest);

for i=1:NumStructs;
    
    dir=[BirdDir DirsOfInterest{i}];
    
    load([dir '/DATSTRUCT.mat']);
    load([dir '/Params.mat']);
   
    typefield=['delay' num2str(Params.HEATMAP.stim_delay)]; % id field
    
    % slide into summary structure
    PARAMSComp.(typefield)=Params; % saves params
    DATSTRUCTComp.(typefield)=DATSTRUCT; 

end
    

%% OLD WAY - manual loading   --------------------------------

% PARAMSComp.(fieldname)=Params; % saves params
% DATSTRUCTComp.(fieldname)=DATSTRUCT; 
% 
% 
% %% SAVES - go to desired summary plot folder.
% save('PARAMSComp','PARAMSComp');
% save('DATSTRUCTComp','DATSTRUCTComp');
% 

% -------------------------------------------------------------

%% PLOT waterfall plot (i.e. all trials, ordered by onset of stim)

delayfields=fieldnames(DATSTRUCTComp);
tPC=PARAMSComp.(delayfields{1}).tf_bins.tPC;


PCcompiled=[];
trigtimes=[];
% 1) compile all trig times and PCs
for i=1:length(delayfields);
    fieldname=delayfields{i};
    
    PCcompiled=[PCcompiled; DATSTRUCTComp.(fieldname).StimNotCatch.PC_ZscoreRelNonstimMean_windowed]; % concatenate all delay types
    trigtimes=[trigtimes PARAMSComp.(fieldname).HEATMAP.TrigTimes_fromEpochStart_delayed];
end


% 2) compile same but, smoothed over trials
PCcompiled_sm=[];
trigtimes_sm=[];
% 1) compile all trig times and PCs
for i=1:length(delayfields);
    fieldname=delayfields{i};
    
    PCcompiled_sm=[PCcompiled_sm; DATSTRUCTComp.(fieldname).StimNotCatch.trialsmoothed.PC_ZscoreRelNonstimMean_windowed]; % concatenate all delay types
    trigtimes_sm=[trigtimes_sm PARAMSComp.(fieldname).HEATMAP.TrigTimes_fromEpochStart_delayed_trialsmoothed'];
end


%% Put things in order

[trigtimes,ind]=sort(trigtimes);
PCcompiled=PCcompiled(ind,:);

[trigtimes_sm,ind]=sort(trigtimes_sm);
PCcompiled_sm=PCcompiled_sm(ind,:);


%% Tack on a number of stim catch trials
% waht fraction of stim trial # to take to sample catch trials

NumStimTrials=size(PCcompiled,1);
NumCatchTrials=round(FracTrialsForCatch*NumStimTrials);

% first cat all catch trials
TMP=[];
for i=1:length(delayfields);
    df=delayfields{i};
    TMP=[TMP; DATSTRUCTComp.(df).StimCatch.PC_ZscoreRelNonstimMean_windowed];
end

% shuffle 
N=size(TMP,1); % num trials
ind=randperm(N);
ind=ind(1:NumCatchTrials);

PCcompiled_catch=TMP(ind,:); % final structs

%% PLOT

% combine catch and nc
Y=[PCcompiled; PCcompiled_catch];

clims=[-3 3];
figure; hold on;
imagesc(tPC,1:size(Y,1),Y,clims)
colormap('spring');
title('Stim trial pitch, z-scored relative to non-stim trial mean and STD, sorted by stim onset');
xlabel('time');
ylabel('trial #')
zlabel('z-score')
colorbar

for i=1:NumStimTrials;
    plot(trigtimes(i)/1000,i,'kx');
end

saveas(gcf,'ZscoredPC_AllTrials_OrderedByStim.fig')
    


NumStimTrials=size(PCcompiled_sm,1);
figure; hold on;
imagesc(tPC,1:NumStimTrials,PCcompiled_sm)
colormap('spring');
title('Stim trial pitch, z-scored relative to non-stim trial mean and STD, sorted by stim onset');
xlabel('time');
ylabel('trial #')
zlabel('z-score')
colorbar

for i=1:NumStimTrials;
    plot(trigtimes_sm(i)/1000,i,'kx');
end

saveas(gcf,'ZscoredPC_AllTrials_OrderedByStim_smoothed.fig')
    


%% NORMALIZE each time bin by itself (i.e. diff)

PCc_mean=mean(PCcompiled,1);

PCcompiled_TimBinDiff=PCcompiled-repmat(PCc_mean,NumStimTrials,1);

clims=[-3 3];
figure; hold on;
imagesc(tPC,1:NumStimTrials,PCcompiled_TimBinDiff)
colormap('spring');
title('Stim trial pitch, z-scored relative to non-stim trial mean and STD, sorted by stim onset');
xlabel('time');
ylabel('trial #')
zlabel('z-score')
colorbar

for i=1:NumStimTrials;
    plot(trigtimes(i)/1000,i,'kx');
end







