%% LT 10/2/15 - PLOT STATS ABOUT THE EXPERIMENTS - E.G. SAMPLE SIZES, DAYS WHEN THINGS HAPPENED, ETC
function [SeqDepPitch_AcrossBirds,Params]= lt_seq_dep_pitch_ACROSSBIRDS_METASTATS(SeqDepPitch_AcrossBirds, Params)

%% PARAMS
NumBirds=length(SeqDepPitch_AcrossBirds.birds);

%% how many birds, experiments, syllables, etc

% ============================================== All experiments
% -------- Initiate figures
count=1;
SubplotsPerFig=4;
subplotrows=2;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];


% --------------- 0) Subplot telling you what class of expeirments are
% being displayed
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
lt_plot_text(0, 0.5, 'ALL EXPERIMENTS');

% -------------- 1) Histogram of number of experiments per bird
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('Number of total experiments'); ylabel('num expts');

X=1:NumBirds;
Y=[];
Xlabels={};

for i=1:NumBirds;
    
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    Xlabels=[Xlabels birdname];
    
    numexpt=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    Y=[Y numexpt];
end

% PLOT
lt_plot_bar(X, Y, {'Color','k'});
set(gca, 'XTick', X);
set(gca, 'XTickLabel', Xlabels);
rotateXLabels(gca, 45)

% TEXT - total experiments
lt_plot_text(1, 2, ['total expts: ' num2str(sum(Y))], 'r');




%% -------------- 2) Histogram of number of syllable (sorted by type) for
% each expt and bird
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('Number of syllables'); ylabel('num syls (red = diff, blue = similar)');

Xlabels={};
SimSyls_all=[];
DiffSyls_all=[];
cnt=0;

for i=1:NumBirds;
    
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        % --- Collect number of similar and different syls
        syls_similar=0;
        syls_different=0;
        
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        for j=1:length(SylsUnique);
            syl=SylsUnique{j};
            
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==0 ...
                    && SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ==1;
                syls_similar=syls_similar+1;
            
            elseif SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).is_target==0 ...
                    && SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ==0;
                syls_different=syls_different+1;
            end
        end
               
        SimSyls_all=[SimSyls_all syls_similar];
        DiffSyls_all=[DiffSyls_all syls_different];
        Xlabels=[Xlabels [birdname(1:4) '-' exptname(end-3:end)]];        
         
        cnt=cnt+1;
        
        % ------ STORE IN STRUCTURE
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.Num_similar_syls=syls_similar;
        SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.Num_diff_syl=syls_different;
        
        
    end
    
    % -- To have a gap between birds
        cnt=cnt+1;
        SimSyls_all=[SimSyls_all nan];
        DiffSyls_all=[DiffSyls_all nan];
        Xlabels=[Xlabels ' '];        
end

% PLOT
X=1:cnt;
lt_plot_bar(X-0.15, SimSyls_all, {'Color','b', 'BarWidth' ,0.3}); % similar
lt_plot_bar(X+0.15, DiffSyls_all, {'Color','r', 'BarWidth' ,0.3}); % diff
set(gca, 'XTick', X);
set(gca, 'XTickLabel', Xlabels);
rotateXLabels(gca, 45)

% TEXT - total syls
lt_plot_text(1, 10, ['total sim syls: ' num2str(nansum(SimSyls_all))], 'b');
lt_plot_text(1, 11, ['total diff syls: ' num2str(nansum(DiffSyls_all))], 'r');



%% ===================== DISPL TARG SYLS AND LEARN DIR FOR EACH BIRD


[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('syls targeted across all expts, by bird');
Birdnames={};
NumUpExpts=0;
NumDnExpts=0;

for i=1:NumBirds;
    
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    StringToPlot=[birdname ':: ' ];
    Birdnames=[Birdnames birdname];
    
    for ii=1:numexperiments
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        targdir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
        switch targdir
            case 1
                targdir='UP';
                NumUpExpts=NumUpExpts+1;
            case -1
                targdir='DN';
                NumDnExpts=NumDnExpts+1;
        end
            
        StringToPlot=[StringToPlot ' | ' targsyl '(' targdir ')'];
    
    end
    
    % ==== plot for this bird
    lt_plot_text(1, i, StringToPlot, 'b');
    
end
ylim([0 NumBirds+1]);
xlim([0 25])

lt_plot_annotation(1, ['up:' num2str(NumUpExpts) '; dn:' num2str(NumDnExpts)], 'k')

%% ------------------------- 3) Experiment durations, timing of certain
% things
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('Number of WN days before start bidir'); ylabel('num WN days');

Xlabel={};
NumWNdays_before_bidir_all=[];
cnt=0;

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        %  ---------- WHEN BEGAN BIDIR LEARNING?
        if ~isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'MultiDir_OneDayBeforeStart_Ind');
            continue;
        end
            
        WNOn_Ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        OneDayBeforeBidir_Ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_OneDayBeforeStart_Ind;
        
        
        % --- skip this expt if bidir was not the first phase (i.e.
        % sometimes samedir was done)
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'SameDir_OneDayBeforeStart_Ind');
            SameDir_OneDayBeforeStart_Ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.SameDir_OneDayBeforeStart_Ind;
            
            if SameDir_OneDayBeforeStart_Ind<OneDayBeforeBidir_Ind
                disp([birdname '-' exptname]);
                continue
                % skip expt
            end
        end
        
        Num_wn_days_before_start_bidir=OneDayBeforeBidir_Ind+1-WNOn_Ind;
        
        
        % --------- OUTPUTS
        Xlabel=[Xlabel [birdname(1:4) '-' exptname(end-3:end)]];
        NumWNdays_before_bidir_all=[NumWNdays_before_bidir_all Num_wn_days_before_start_bidir];
        
        cnt=cnt+1;
        
    end
end

% PLOT
X=1:cnt;

lt_plot_bar(X, NumWNdays_before_bidir_all);

set(gca, 'XTick', X);
set(gca, 'XTickLabel', Xlabel);
rotateXLabels(gca, 45)

lt_plot_annotation(1, ['mean=' num2str(mean(NumWNdays_before_bidir_all)) '; std=' num2str(std(NumWNdays_before_bidir_all))]);


%% ==== distribution of num days during which increased threshold for single-context training

[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('Number of days adapting threshold'); 
xlabel('last WN day with thresh update (starting count from good days');
NumDaysOK=[]; % last day substantial change
NumDaysComplete=[]; % last day for any change, even small.


for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        % --- number of days
        numdays=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.NumDaysAdapticeThr_col1OK_col2complete;
        
        NumDaysOK=[NumDaysOK numdays(1)];
        NumDaysComplete=[NumDaysComplete numdays(2)];
        
        
    end
end


lt_plot_cdf(NumDaysOK, 'b', 1);
lt_plot_cdf(NumDaysComplete, 'k', 1);
lt_plot_annotation(1, ['blue=OK;  black=complete'], 'k');



%% ====== number of motifs per bird, distribution of num syls in motif
NumMotifsPerBird=[];
NumMotifsPerBird_regexp=[];

NumSylsPerMotif=[];
NumSylsPerMotif_regexp=[];
for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexperiments=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        % === number of motifs
        nummotifs=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.FieldsInOrder);
        nummotifs_regexp=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.RegularExpressionsList);
        
        NumMotifsPerBird=[NumMotifsPerBird nummotifs];
        NumMotifsPerBird_regexp=[NumMotifsPerBird_regexp nummotifs_regexp];
        
        % === number of syllables in motif
        for j=1:nummotifs
            numsyls=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.FieldsInOrder{j});
                      NumSylsPerMotif=[NumSylsPerMotif numsyls ];  
        end
        
        for j=1:nummotifs_regexp
            numsyls2=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.RegularExpressionsList{j});
                      NumSylsPerMotif_regexp=[NumSylsPerMotif_regexp numsyls2 ];  
        end
    end
end

[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('num motifs per bird (rd=syllists); (bk=regexp)');

lt_plot_cdf(NumMotifsPerBird, 'r', 1)
lt_plot_cdf(NumMotifsPerBird_regexp, 'k', 1)

lt_plot_annotation(1, ['mean, std=' num2str(mean(NumMotifsPerBird_regexp)) ', ' num2str(std(NumMotifsPerBird_regexp))], 'k')
 lt_plot_annotation(2, ['total = ' num2str(sum(NumMotifsPerBird_regexp))], 'k')
   
    
[fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
title('num syls per motif (rd=syllist); bk=regexp');

lt_plot_cdf(NumSylsPerMotif, 'r', 1)
lt_plot_cdf(NumSylsPerMotif_regexp, 'k', 1)

lt_plot_annotation(1, ['mean, std=' num2str(mean(NumSylsPerMotif_regexp)) ', ' num2str(std(NumSylsPerMotif_regexp))], 'k')
    




%% ============================================== Expeirments with single
% dir



% =========================================== Bidir experiments




% ============================================= LMAN only


