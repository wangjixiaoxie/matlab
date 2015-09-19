clear all; close all;
% ENTER PATHS TO FIND THE COMPILED DATA
% have to run lt_compile_seq_dep_pitch_data_PLOTDirUndir first for each bird.  That saves a data structure. Enter the path of that structure below.

% experiment where started with one, then hit two.
% WILL ONLY LOOK AT LEARNING DURING ONE TARGET
SeqDepPitch_AcrossBirds.birds.pu11wh87_3.path=...
    '/bluejay3/lucas/birds/pu11wh87/compile_seq_dep_pitch_data_SeqDepPitchShift2_PLOTS/18Nov2014to20Dec2014_made06Jan2015_1410/SeqDepPitch.mat';

% SeqDepPitch_AcrossBirds.birds.pu37wh20_2.path=...
%     '/bluejay3/lucas/birds/pu37wh20/compile_seq_dep_pitch_data_SeqDepPitchShift_PLOTS/26Oct2014to18Nov2014_made03Dec2014_2146/SeqDepPitch';

SeqDepPitch_AcrossBirds.birds.pu37wh20_1.path=...
    '/bluejay3/lucas/birds/pu37wh20/compile_seq_dep_pitch_data_SeqDepPitchShift_PLOTS/26Oct2014to22Nov2014_made14Jan2015_2225/SeqDepPitch.mat';

SeqDepPitch_AcrossBirds.birds.pu53wh88_3.path=...
        '/bluejay3/lucas/birds/pu53wh88/compile_seq_dep_pitch_data_SeqDepPitchShift_PLOTS/26Oct2014to18Nov2014_made11Jan2015_1611/SeqDepPitch.mat';
    
if (0); % ignore pu64 because he only has similar syls. Is useful when looking at distance effects. Also, has 2 learning epochs.
    SeqDepPitch_AcrossBirds.birds.pu64bk13_1.path=...
        '/bluejay3/lucas/birds/pu64bk13/compile_seq_dep_pitch_data_SeqDepPitchShift_PLOTS/10Dec2014to28Dec2014_made05Jan2015_2248/SeqDepPitch.mat';
end


%% Compile into one structure

BirdList=fieldnames(SeqDepPitch_AcrossBirds.birds);
NumBirds=length(BirdList);

for i=1:length(BirdList);
    bird=BirdList{i};
    SeqDepPitch_AcrossBirds.birds.(bird)=load(SeqDepPitch_AcrossBirds.birds.(bird).path);
end


% MAKE edits:

try
    % remove bccB as this was actually also a target (along with dccB);
    ii=find(strcmp(SeqDepPitch_AcrossBirds.birds.pu37wh20_1.FinalStruct.params.ARGUMENTS.SylLists.SylsSame,'bccB'));
    SeqDepPitch_AcrossBirds.birds.pu37wh20_1.FinalStruct.params.ARGUMENTS.SylLists.SylsSame(ii)=[];
    ii=find(strcmp(SeqDepPitch_AcrossBirds.birds.pu37wh20_1.FinalStruct.params.ARGUMENTS.SylLists.SylsSame,'qccB'));
    SeqDepPitch_AcrossBirds.birds.pu37wh20_1.FinalStruct.params.ARGUMENTS.SylLists.SylsSame(ii)=[];
catch err
    disp('pu37 edit issue- look');
end



%% PARAMS

Params.PlotColors=lt_make_plot_colors(NumBirds,0);



%% Check direction of learning of target

for kk=1:NumBirds;
    birdfield=BirdList{kk};
    targsyl=SeqDepPitch_AcrossBirds.birds.(birdfield).FinalStruct.params.ARGUMENTS.SylLists.TargetSyls{1};
    
    DirOfShift.(birdfield)=sign(SeqDepPitch_AcrossBirds.birds.(birdfield).FinalStruct.EpochData.WNdaysSlidingWin{end}.UNDIR.(targsyl).meanFF_minusBaseline);
    
    % make sure learning greater than 1.5 z-score, or else query user
    if abs(SeqDepPitch_AcrossBirds.birds.(birdfield).FinalStruct.EpochData.WNdaysSlidingWin{end}.UNDIR.(targsyl).Zscore)<1.5;
        disp(['Potential problem - last day of learning for ' birdfield ' low, and useing that day to determine direction of learning.']);
    end
end

disp('Direction of pitch shift:');
disp(DirOfShift)



% 2) Normalize learning by maximum learning of target syl


%% PLOT PEAK OF LEARNING - e.g. 3 day mean
targsyl=[];
for i=1:NumBirds;
    birdfield=BirdList{i};
    
    % learning at targ syl
    targsyl=SeqDepPitch_AcrossBirds.birds.(birdfield).FinalStruct.params.ARGUMENTS.SylLists.TargetSyls{:};
    % final learning vals
    if strcmp(birdfield,'pu11wh87_3')==1; % only use start when driving one syl, not two.
        FinalLearningVals.(birdfield)=SeqDepPitch_AcrossBirds.birds.pu11wh87_3.FinalStruct.EpochData.Snapshot.from09Dec2014to11Dec2014.UNDIR;
    else % all other cases
        FinalLearningVals.(birdfield)=SeqDepPitch_AcrossBirds.birds.(birdfield).FinalStruct.EpochData.WNdaysSlidingWin{end}.UNDIR; % structure with final learning (final 3 days of WN);
    end
    
    % Normalize all
    NormFactor=FinalLearningVals.(birdfield).(targsyl).meanFF_minusBaseline; % learning at target
    
    % divide all by target leanring
    SylFields=fieldnames(FinalLearningVals.(birdfield));
    for j=1:length(SylFields);
        sylfield=SylFields{j};
        FinalLearningVals.(birdfield).(sylfield).NormalizedToTarg_hz=...
            FinalLearningVals.(birdfield).(sylfield).meanFF_minusBaseline/NormFactor;
    end
end

FinalLearningVals;


% PLOT
figure; hold on;
Ymean=[];
Ystd=[];
Ysem=[];


for i=1:NumBirds;
    birdfield=BirdList{i};
    
    % what are sim and diff syls?
    SimSyls=SeqDepPitch_AcrossBirds.birds.(birdfield).FinalStruct.params.ARGUMENTS.SylLists.SylsSame;
    DiffSyls=SeqDepPitch_AcrossBirds.birds.(birdfield).FinalStruct.params.ARGUMENTS.SylLists.SylsDifferent;
    TargSyl=SeqDepPitch_AcrossBirds.birds.(birdfield).FinalStruct.params.ARGUMENTS.SylLists.TargetSyls{:};
    
    % PLOT SIMILAR
    subplot(1,3,1); hold on;
    xlim([0.8-0.02*NumBirds 1.3+0.02*NumBirds]);
    set(gca,'XTick',[0.85 1.0 1.2],'XTickLabel',{'All Birds','All Syls','Syl Name'});
    title('Similar syls (mean/SEM)');
    
    Y=[];
    for j=1:length(SimSyls);
        sylfield=SimSyls{j};
        Y=[Y FinalLearningVals.(birdfield).(sylfield).NormalizedToTarg_hz]; % compile datapoints
    end
    plot(1+0.02*i,Y,'.','Color',Params.PlotColors{i});
    text(ones(length(Y),1)+0.07+(i+NumBirds)*0.02,Y,SimSyls,'Color',Params.PlotColors{i})
    
    % plot targ syl different color
    %     plot(1+0.02*i,FinalLearningVals.(birdfield).(TargSyl).NormalizedToTarg_hz,'o','Color',Params.PlotColors{i},...
    %         'MarkerFaceColor',Params.PlotColors{i});
    % plot means for each bird
    Ymean.sim(i)=mean(Y);
    Ystd.sim(i)=std(Y);
    Ysem.sim(i)=Ystd.sim(i)/sqrt(length(Y)-1);
    h1(i)=errorbar(0.9-0.02*i,Ymean.sim(i),Ysem.sim(i),'d','Color',Params.PlotColors{i},'MarkerFaceColor',Params.PlotColors{i},'MarkerSize',7);
    
    % save individual syl values
    Yvals.sim.(birdfield)=Y;
    
    % line at y=0;
    line(xlim,[0 0],'LineStyle','--','Color','k')
    
    
    
    % PLOT DIFFERENT
    subplot(1,3,2); hold on;
    xlim([0.8-0.02*NumBirds 1.3+0.02*NumBirds]);
    title('Different syls (mean/SEM)');
    set(gca,'XTick',[0.85 1.0 1.2],'XTickLabel',{'All Birds','All Syls','Syl Name'});
    
    Y=[];
    for j=1:length(DiffSyls);
        sylfield=DiffSyls{j};
        Y=[Y FinalLearningVals.(birdfield).(sylfield).NormalizedToTarg_hz]; % compile datapoints
    end
    plot(1+0.02*i,Y,'.','Color',Params.PlotColors{i});
    text(ones(length(Y),1)+0.1+NumBirds*0.02,Y,DiffSyls,'Color',Params.PlotColors{i})
    
    % plot means for each bird
    Ymean.dif(i)=mean(Y);
    Ystd.dif(i)=std(Y);
    Ysem.dif(i)=Ystd.dif(i)/sqrt(length(Y)-1);
    h1(i)=errorbar(0.9-0.02*i,Ymean.dif(i),Ysem.dif(i),'d','Color',Params.PlotColors{i},'MarkerFaceColor',Params.PlotColors{i},'MarkerSize',7);
    
    % save individual syl values
    Yvals.dif.(birdfield)=Y;
    % line at y=0;
    line(xlim,[0 0],'LineStyle','--','Color','k')
    
    
    
    
    % THIRD, Plot on same plot - same y scale
    subplot(1,3,3); hold on;
    xlim([0.2 2.8]);
    title('Similar and Different syls combined (with All Syls mean/SEM)');
    set(gca,'XTick',[1 2],'XTickLabel',{'Similar Syls','Diff Syls'});
    
    
    % Similar
    plot(1+0.02*i,Yvals.sim.(birdfield),'.','Color',Params.PlotColors{i});
    
    Y=Yvals.sim.(birdfield);
    plot(1+NumBirds*0.02+0.02*(i+2),Ymean.sim(i),'d','Color',Params.PlotColors{i},'MarkerFaceColor',Params.PlotColors{i}),'MarkerSize',7;

    
    % Different
    plot(2+0.02*i,Yvals.dif.(birdfield),'.','Color',Params.PlotColors{i});
    plot(2+NumBirds*0.02+0.02*(i+2),Ymean.dif(i),'d','Color',Params.PlotColors{i},'MarkerFaceColor',Params.PlotColors{i},'MarkerSize',7);

    line(xlim,[0 0],'LineStyle','--','Color','k')
    
    
end


legend(h1,BirdList);

% Plot mean across all birds - FOR THE SIMILAR AND DIFF ALONE PLOTS
% SIMILAR - using birds
YYmean=mean(Ymean.sim);
YYstd=std(Ymean.sim);
YYsem=YYstd/sqrt(length(Ymean.sim));

subplot(1,3,1); hold on;
errorbar(0.9-0.02*NumBirds-0.04,YYmean,YYsem,'s','Color','k','MarkerSize',8);
% SIMILAR - using indivudal syls
YY=[];
for i=1:NumBirds;
    birdfield=BirdList{i};
    YY=[YY Yvals.sim.(birdfield)];
end
YYmean=mean(YY);
YYstd=std(YY);
YYsem=YYstd/sqrt(length(YY));

subplot(1,3,1); hold on;
errorbar(1,YYmean,YYsem,'s','Color','k','MarkerSize',8);
% for combined plot
subplot(1,3,3);
errorbar(1,YYmean,YYsem,'s','Color','k','MarkerSize',8);


% DIFFERENT - using birds
YYmean=mean(Ymean.dif);
YYstd=std(Ymean.dif);
YYsem=YYstd/sqrt(length(Ymean.dif));

subplot(1,3,2); hold on;
errorbar(0.9-0.02*NumBirds-0.04,YYmean,YYsem,'s','Color','k','MarkerSize',8);
% DIFFERENT - using indivudal syls
YY=[];
for i=1:NumBirds;
    birdfield=BirdList{i};
    YY=[YY Yvals.dif.(birdfield)];
end
YYmean=mean(YY);
YYstd=std(YY);
YYsem=YYstd/sqrt(length(YY));

subplot(1,3,2); hold on;
errorbar(1,YYmean,YYsem,'s','Color','k','MarkerSize',8);
% for combined plot
subplot(1,3,3);
errorbar(2,YYmean,YYsem,'s','Color','k','MarkerSize',8);


%% PLOT DAY BY DAY DATA

% PREPROCESSING

% 1) Number of overlapping WN days across birds
XX=[]; % num of days having WN
for i=1:NumBirds;
    birdfield=BirdList{i};
    
    % PLACE EXCEPTIONS HERE - i.e. last DAY of WN as listed is not accurate
    %, since might have driven one syl, then another (or 2) in same
    % experiment
    
    if strcmp(birdfield,'pu11wh87_3')==1; % first drove one, then did bidirectional.
        DaysWithWN.(birdfield).FirstWNDay=...
            SeqDepPitch_AcrossBirds.birds.pu11wh87_3.FinalStruct.params.WNTimeOnInd; % 1st day of WN
        DaysWithWN.(birdfield).FirstNoWNDay=25; % 1st day of no WN
        
        DaysWithWN.(birdfield).NumDaysWithWN=DaysWithWN.(birdfield).FirstNoWNDay-DaysWithWN.(birdfield).FirstWNDay;
        
        XX=[XX, DaysWithWN.(birdfield).NumDaysWithWN]; % to figure out minimum number of days across birds;
        
    else
        DaysWithWN.(birdfield).FirstWNDay=...
            SeqDepPitch_AcrossBirds.birds.(birdfield).FinalStruct.params.WNTimeOnInd; % 1st day of WN
        DaysWithWN.(birdfield).FirstNoWNDay=...
            SeqDepPitch_AcrossBirds.birds.(birdfield).FinalStruct.params.WNTimeOffInd; % 1st day of no WN
        
        DaysWithWN.(birdfield).NumDaysWithWN=DaysWithWN.(birdfield).FirstNoWNDay-DaysWithWN.(birdfield).FirstWNDay;
        
        XX=[XX, DaysWithWN.(birdfield).NumDaysWithWN]; % to figure out minimum number of days across birds;
    end
end


MaxOverlapWNDays=min(XX);
DaysWithWN;



% 2) EXTRACT DAY BY DAY DATA
% for i=1:NumBirds;
%     birdfield=BirdList{i};
%     
%     
% %     DaysToPlot.StartLearning.(birdfield)=[DaysWithWN.pu11wh87_3.NumDaysWithWN
% 
% 
% % 1) Plot day by day average - same syl and diff. syl
% 
% 
% 
% SeqDepPitch_AcrossBirds.birds.pu37wh20_1.FinalStruct.DataMatrix.MatrixOverDaysforSylLists.SylsSame.UNDIR.meanFF_minusBaseline


%% PLOT - start of learning, across maximum num of overlap days between
% birds.
figure; hold on;
Y=[];

for i=1:NumBirds;
birdfield=BirdList{i};
x1=DaysWithWN.(birdfield).FirstWNDay-3; % first wn day, minus 3 to include 3 baseline days.
x2=(x1+3)+MaxOverlapWNDays-1; % last data day

% TARGSYL
TargSyl=SeqDepPitch_AcrossBirds.birds.(birdfield).FinalStruct.params.ARGUMENTS.SylLists.TargetSyls{:};
TargDataMatrix=SeqDepPitch_AcrossBirds.birds.(birdfield).FinalStruct.DataMatrix.UNDIR.(TargSyl).meanFF_DevFromBase; % all days leanring
TargDataMatrix=TargDataMatrix(x1:x2); % only 3 days baseline + WN days

% invert data if is driving pitch down
TargDataMatrix=TargDataMatrix.*DirOfShift.(birdfield);

% normalization factor
NormFactor=max(max(TargDataMatrix));
TargDataMatrix_norm=TargDataMatrix./NormFactor;


% SIMILAR SYLABLES
subplot(1,2,1); hold on;
title('Similar Syllables - syl as datapoint');
ylabel('Learning, normalized to max learning at target');
xlabel('Days');

% extract data
Y=SeqDepPitch_AcrossBirds.birds.(birdfield).FinalStruct.DataMatrix.MatrixOverDaysforSylLists.SylsSame.UNDIR.meanFF_minusBaseline(x1:x2,:);

% invert data if is driving pitch down
Y=Y.*DirOfShift.(birdfield);
% normalize each day by maximum learning
Y=Y./NormFactor;

% PLOT
plot(Y,'-','Color',Params.PlotColors{i});
line([3.5 3.5],ylim,'Color','k','LineStyle','--');
% put into structure to get mean
Yvals_daybyday.sim.(birdfield)=Y;
% plot targ syl
plot(TargDataMatrix_norm,'--','Color',Params.PlotColors{i});



% DIFFERENT SYLLABLES
Y=[];
subplot(1,2,2); hold on;
title('Different Syllables - syl as datapoint');
ylabel('Learning, normalized to max learning at target');
xlabel('Days');

% extract data
Y=SeqDepPitch_AcrossBirds.birds.(birdfield).FinalStruct.DataMatrix.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.meanFF_minusBaseline(x1:x2,:);

% invert data if is driving pitch down
Y=Y.*DirOfShift.(birdfield);
% normalize each day by maximum learning
Y=Y./NormFactor;

% PLOT
plot(Y,'-','Color',Params.PlotColors{i});
line([3.5 3.5],ylim,'Color','k','LineStyle','--');
% put into structure to get mean
Yvals_daybyday.dif.(birdfield)=Y;
% plot targ syl
plot(TargDataMatrix_norm,'--','Color',Params.PlotColors{i});
end

% GET MEAN (for each day)
% Similar syls
subplot(1,2,1); hold on;
YY=[];
for i=1:NumBirds;
    birdfield=BirdList{i};
    YY=[YY Yvals_daybyday.sim.(birdfield)]; % rows days, cols syls
end
Ymean.sim=nanmean(YY,2);
Ystd.sim=nanstd(YY,0,2);
for i=1:size(YY,1); % get sample size - diff days might have diff N
    N(i)=length(YY(i,:))-sum(isnan(YY(i,:))); % total number of syls, minus number that have NAN.
end
Ysem.sim=Ystd.sim./sqrt(N-1);
% Plot
errorbar(Ymean.sim,Ysem.sim,'-o','Color','k')

% Diff Syls
subplot(1,2,2); hold on;
YY=[];
for i=1:NumBirds;
    birdfield=BirdList{i};
    YY=[YY Yvals_daybyday.dif.(birdfield)]; % rows days, cols syls
end
Ymean.dif=nanmean(YY,2);
Ystd.dif=nanstd(YY,0,2);
for i=1:size(YY,1); % get sample size - diff days might have diff N
    N(i)=length(YY(i,:))-sum(isnan(YY(i,:))); % total number of syls, minus number that have NAN.
end
Ysem.dif=Ystd.dif./sqrt(N-1);
% Plot
errorbar(Ymean.dif,Ysem.dif,'-o','Color','k')



% PLOT with each day as bird - mean and std
figure; hold on; 
Ybird=struct;

% Similar syls
subplot(1,2,1); hold on;
title('Similar Syls - bird as datapoint)');
ylabel('Learning - normalized target max (black: mean/SEM)');
xlabel('Days');

for i=1:NumBirds;
    birdfield=BirdList{i};
    Ybird.sim.vals(:,i)=mean(Yvals_daybyday.sim.(birdfield),2); % each bird is one col, rows are days
end

% get sample size for each day
N=sum(~isnan(Ybird.sim.vals),2);  
% get stats
Ybird.sim.mean=nanmean(Ybird.sim.vals,2);
Ybird.sim.std=nanstd(Ybird.sim.vals,0,2);
Ybird.sim.sem=Ybird.sim.std./sqrt(N-1);
% plot each bird
for i=1:NumBirds;
    h2(i)=plot(Ybird.sim.vals(:,i),'Color',Params.PlotColors{i});
end
legend(h2,BirdList);
% plot stats across birds
errorbar(Ybird.sim.mean,Ybird.sim.sem,'o-k');
% line for baseline
line([3.5 3.5],ylim,'Color','k','LineStyle','--');


% Different syls
subplot(1,2,2); hold on;
title('Different Syls - bird as datapoint)');
xlabel('Days');

for i=1:NumBirds;
    birdfield=BirdList{i};
    Ybird.dif.vals(:,i)=mean(Yvals_daybyday.dif.(birdfield),2); % each bird is one col, rows are days
end

% get sample size for each day
N=sum(~isnan(Ybird.dif.vals),2);  
% get stats
Ybird.dif.mean=nanmean(Ybird.dif.vals,2);
Ybird.dif.std=nanstd(Ybird.dif.vals,0,2);
Ybird.dif.sem=Ybird.dif.std./sqrt(N-1);
% plot each bird
for i=1:NumBirds;
    h2(i)=plot(Ybird.dif.vals(:,i),'Color',Params.PlotColors{i});
end
legend(h2,BirdList);
% plot stats across birds
errorbar(Ybird.dif.mean,Ybird.dif.sem,'o-k');
% line for baseline
line([3.5 3.5],ylim,'Color','k','LineStyle','--');

    
% PLOT JUST MEAN
figure; hold on;
% Syllable as data point
% --Similar syls
hplot1(1)=subplot(1,2,1); hold on;
title('Syllable as datapoint');
ylabel('Mean learning (normalized to targ) +/- SEM');
xlabel('Days');
h3(1)=errorbar(Ymean.sim,Ysem.sim,'-o','Color','b');
% --Diff Syls
h3(2)=errorbar(Ymean.dif,Ysem.dif,'-o','Color','r');
line([3.5 3.5],ylim,'Color','k','LineStyle','--');
line(xlim,[0 0],'Color','k','LineStyle','--');

legend(h3,{'Similar','Different'});

% Bird as datapoint
% --Similar syls
hplot1(2)=subplot(1,2,2); hold on;
title('Bird as datapoint');
xlabel('Days');
errorbar(Ybird.sim.mean,Ybird.sim.sem,'o-b');
% --Diff Syls
errorbar(Ybird.dif.mean,Ybird.dif.sem,'o-r');
line([3.5 3.5],ylim,'Color','k','LineStyle','--');
line(xlim,[0 0],'Color','k','LineStyle','--');
linkaxes(hplot1,'xy')
    


% PLOT VARIABILITY - 

% % PLOT VARIABILITY - ACROSS SYLS
% % Similar syls
% % -- get mean FF to quantify CV
% for i=1:NumBirds;
%     birdfield=BirdList{i};
%     
%     
% 
% end
% 
%     
%     Yvals_daybyday.sim.(birdfield)
% Ymean_abs=
% Ycv=Ystd.sim./Ymean.sim
% 




%% PLOT OVER CONSOLIDATION DAYS






%% PLOT and Compile data.

% SAME SYLS
FieldToCheck='SylsSame';

% PREP WORK
% 1) Find the birds that have this field

XX=[];
for kk=1:NumBirds;
    birdfield=BirdList{kk};
    if ~isfield(SeqDepPitch_AcrossBirds.birds.(birdfield).FinalStruct.params.ARGUMENTS.SylLists,FieldToCheck); % then this bird is thrown out
        XX=[XX, kk];
    end
end

BirdListNew=BirdList;
BirdListNew(XX)=[];
NumBirdsNew=length(BirdListNew);



% PLOT
figure; hold on;
subplot(1,2,1); hold on;

Yall=[];
birdcolors=lt_make_plot_colors(NumBirdsNew,0);
for kk=1:NumBirdsNew;
    
    birdfield=BirdListNew{kk};
    
    % PLOT SIMILAR SYLS - early and late
    FieldsList=SeqDepPitch_AcrossBirds.birds.(birdfield).FinalStruct.params.ARGUMENTS.SylLists.SylsSame;
    
    for ii=1:length(FieldsList); % how many fields within this set?
        syl=FieldsList{ii}; % actual syl name (e.g. 'a')
        
        Y=[SeqDepPitch_AcrossBirds.birds.(birdfield).FinalStruct.EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.meanFF_minusBaseline(1,ii),...
            SeqDepPitch_AcrossBirds.birds.(birdfield).FinalStruct.EpochData.MatrixOverDaysforSylLists.SylsSame.UNDIR.meanFF_minusBaseline(end,ii)];
        
        % make all shifts go in positive direction
        Y=Y.*DirOfShift.(birdfield);
        
        plot(Y,'-ob');
        
        Xlims=xlim;
        text(Xlims(2)+0.1,Y(2),syl,'Color',birdcolors{kk}); % label the datapoint
        
        % compile, to get mean across all syls and birds.
        Yall=[Yall;Y];
        
    end
end

% Calculate mean
Ymean=mean(Yall,1);
Ysd=std(Yall,0,1);
Ysem=Ysd/sqrt(size(Yall,1)-1);

% plot mean and SEM
X=[1 2];
Xtl={'Early','Late'};
errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','b');

% Format plot
title('Similar Syllables');
ylabel('Pitch (minus baseline) (hz)');
set(gca,'XTick',X);
set(gca,'XTickLabel',Xtl);

Ylimits=ylim;
xlim([0.5 2.5])
lt_plot_zeroline;



% PLOT DIFF SYLS
subplot(1,2,2); hold on;
Yall=[];
for kk=1:NumBirds;
    birdfield=BirdList{kk};
    try
        % PLOT DIFFERENT SYLS
        FieldsList=SeqDepPitch_AcrossBirds.birds.(birdfield).FinalStruct.params.ARGUMENTS.SylLists.SylsDifferent;
        
        for ii=1:length(FieldsList); % how many fields within this set?
            syl=FieldsList{ii}; % actual syl name (e.g. 'a')
            
            Y=[SeqDepPitch_AcrossBirds.birds.(birdfield).FinalStruct.EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.meanFF_minusBaseline(1,ii),...
                SeqDepPitch_AcrossBirds.birds.(birdfield).FinalStruct.EpochData.MatrixOverDaysforSylLists.SylsDifferent.UNDIR.meanFF_minusBaseline(end,ii)];
            
            % make all shifts go in positive direction
            Y=Y.*DirOfShift.(birdfield);
            
            plot(Y,'-or');
            
            Xlims=xlim;
            text(Xlims(2)+0.1,Y(2),syl,'Color',birdcolors{kk}); % label the datapoint
            
            % compile, to get mean across all syls and birds.
            Yall=[Yall;Y];
            
        end
        
    catch err
        disp([birdfield ' does not have SameSyl fields - skipping him']);
        continue
    end
end

% Calculate mean
Ymean=mean(Yall,1);
Ysd=std(Yall,0,1);
Ysem=Ysd/sqrt(size(Yall,1)-1);

% plot mean and SEM
X=[1 2];
Xtl={'Early','Late'};
errorbar(X-0.1,Ymean,Ysem,'-ok','MarkerFaceColor','r');

% Format plot
title('Different Syllables');
ylabel('Pitch (minus baseline) (hz)');
set(gca,'XTick',X);
set(gca,'XTickLabel',Xtl);
lt_plot_zeroline;

ylim(Ylimits);
xlim([0.5 2.5])
lt_plot_zeroline;

% GLOBAL title
subtitle('Mean learning (3-day bins) - All syls across all birds');


%% PLOT DAY BY DAY
% will align all birds by start, and plto means for all days that have
% data from all birds.



