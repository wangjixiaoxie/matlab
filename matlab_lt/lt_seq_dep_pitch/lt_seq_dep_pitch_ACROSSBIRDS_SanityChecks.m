%% LT script written 11/25/15 - 

%% MAKE SURE BIDIR START AND END INDS ARE IDENTICAL FOR PLOT LEARNING AND THIS CODE

for i=1:length(SeqDepPitch_AcrossBirds.birds)
    
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    for ii=1:numexpts;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        disp(' ---');
        disp([birdname '-' exptname]);
        
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning, 'timeline');
            
            if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline, 'bidir_start_Inds');
                
                
                % ==== does it have bidir
                BidirBeforeStart1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_OneDayBeforeStart_Ind;
                BidirBeforeStart2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.timeline.bidir_start_Inds;
                
                disp(['start: ' num2str(BidirBeforeStart1+1) '--' num2str(BidirBeforeStart2)]);
                
            end
        end
    end
end


    


    
    

%% CHECK THAT SEQUENCE STATS ARE CORRECT - E.G. PRESYL, POSTSYL, ETC

NumBirds=length(SeqDepPitch_AcrossBirds.birds);
disp(' ');
disp('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++');


for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexperiments =length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        disp(' ');
        disp(['===========' birdname '-' exptname]);
        
targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
presyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).preceding_syl;
twosylback=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).two_syl_back;
singlesyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targsyl).single_syl;

if isnan(presyl);
    presyl= 'NA';
end
if isnan(twosylback);
    twosylback='NA';
end

disp(['TARGET: ' twosylback ' - ' presyl ' - ' singlesyl ' == ' targsyl ]);

        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            if strcmp(syl, targsyl)
                continue
            end
            
          
            
% ==== extract presyl, postsyl, and whether is similar to targ

presyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).preceding_syl;
twosylback=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).two_syl_back;
singlesyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).single_syl;
   
if isnan(presyl);
    presyl= 'NA';
end
if isnan(twosylback);
    twosylback='NA';
end

sim=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).similar_to_targ;
presim=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).presyl_similar_to_targ_presyl;
twopresim=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).two_syl_back_same_as_targ;


disp(['SYL: ' twosylback ' - ' presyl ' - ' singlesyl ' == ' syl ...
    ' == sim (' num2str(sim) '); presim (' num2str(presim) '); twopresim (' num2str(twopresim) ')']);

        end
    end
end

%% === compare reg exp motifs with fields motifs - ideally they are the same

NumBirds=length(SeqDepPitch_AcrossBirds.birds);

for i=1:NumBirds
    numexpts = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
   
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    for ii=1:numexpts
   exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
   disp(' ');
   disp([' === ' birdname ' -- ' exptname]);
   
        % ==== reg expt motifs
        RegExps=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.Params.RegExpr.expressions;
        
        FieldsInOrder=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.SylLists.FieldsInOrder;
                
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.SylFields_Unique;
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;

        target=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        
        ind1=find(strcmp(PARAMS.global.SylsToRemove, birdname));
        ind2=find(strcmp(PARAMS.global.SylsToRemove, exptname));
        ind3=intersect(ind1+1, ind2);
        
        % ==== display
        disp(['regexp: ' RegExps ]);
        disp(['syl fields: ' FieldsInOrder{:}]);
        disp(['target: ' {target}]);
        disp(['unique syls: ' SylsUnique]);
       
        % display syls removed if they exist
        if ~isempty(ind3)
            disp(['(syls removed): ' PARAMS.global.SylsToRemove{ind3+1}]);
        end
        
        
        
        
    end
end

%% NOT ENOUGH SAMPLE FOR ACOUSTIC? 
% Display sample sizes for all expts and syls for structure stats (passing
% duration threshold), and all plot learning (all)
NumBirds=length(SeqDepPitch_AcrossBirds.birds);
disp(' ');
disp(' ======================================= ');
disp('N_acoustic/N_all [NOT MUSC DATA]');
disp('------------');

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexperiments =length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.SylFields_Unique;
        
        disp(' ');
        disp([birdname '-' exptname]);
        
        
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
            N_acoustic=size(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_all_rends,1);
            N_all=numel(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).rawFF_WithinTimeWindow);
                
            else
                
            N_acoustic=size(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data.(syl).fv_baseline_all_rends,1);
            N_all=numel(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).rawFF);
            
            end
            
            if N_acoustic<N_all;
                disp([syl ':  ' num2str(N_acoustic) '/' num2str(N_all)]);
            else
                disp(['-- GOOD -- ' syl ':  ' num2str(N_acoustic) '/' num2str(N_all)]);
            end
            
            
            
        end
    end
end
        

%% NOT ENOUGH SAMPLE FOR ACOUSTIC? [LMAN DATA]
% Display sample sizes for all expts and syls for structure stats (passing
% duration threshold), and all plot learning (all)
NumBirds=length(SeqDepPitch_AcrossBirds.birds);
disp(' ');
disp(' ======================================= ');
disp('N_acoustic/N_all [MUSC data]');
disp('------------');

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexperiments =length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==0
            continue
        end
        
        
        syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.SylFields_Unique;
        
        disp(' ');
        disp([birdname '-' exptname]);
        
        
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            N_acoustic=size(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_StructureStats.data_MUSC.(syl).fv_baseline_all_rends,1);
            N_all=numel(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).rawFF_WithinTimeWindow);
                
            
            if N_acoustic<N_all;
                disp([syl ':  ' num2str(N_acoustic) '/' num2str(N_all)]);
            else
                disp(['-- GOOD -- ' syl ':  ' num2str(N_acoustic) '/' num2str(N_all)]);
            end
            
            
            
        end
    end
end
        

%% NOT ENOUGH SAMPLE FOR CORRELATIONS? [REGEXP] [PBS DATA]
NumBirds=length(SeqDepPitch_AcrossBirds.birds);
disp(' ');
disp(' ======================================= ');
disp('N_regexp/N_all [NOT MUSC DATA]');
disp('------------');

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexperiments =length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        disp(' ');
        disp([birdname '-' exptname]);
        
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            % -- get regexpt motif
            motifnum=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).regexp_motifnum;
            posinmotif=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).regexp_PosInMotif_thissyl;
            
            if isempty(motifnum)
                disp([syl ': EMPTY (not in any regexp motif']);
                continue
            end
            
            N_regexp=sum(~isnan(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline.data_WithOutlier{motifnum}.FFvals(:, posinmotif))); % all samples for this syl that are not nan
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1
                N_all=numel(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).rawFF_WithinTimeWindow);
            else
                N_all=numel(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData.Baseline.(syl).rawFF);
            end
            
            if N_regexp<N_all;
                disp(['[motif' num2str(motifnum) '|pos' num2str(posinmotif) '] ' syl ': ' num2str(N_regexp) '/' num2str(N_all)]);
            else
                disp(['[motif' num2str(motifnum) '|pos' num2str(posinmotif) '] -- GOOD -- ' syl ':  ' num2str(N_regexp) '/' num2str(N_all)]);
            end
            
            
            
        end
    end
end
        

%% NOT ENOUGH SAMPLE FOR CORRELATIONS? [REGEXP] [MUSC DATA]
NumBirds=length(SeqDepPitch_AcrossBirds.birds);
disp(' ');
disp(' ======================================= ');
disp('N_regexp/N_all [NOT MUSC DATA]');
disp('------------');

for i=1:NumBirds;
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexperiments =length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexperiments;
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==0
            continue
        end
        
        syls_unique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        
        disp(' ');
        disp([birdname '-' exptname]);
        
        for j=1:length(syls_unique);
            syl=syls_unique{j};
            
            % -- get regexpt motif
            motifnum=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).regexp_motifnum;
            posinmotif=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(syl).regexp_PosInMotif_thissyl;
            
            if isempty(motifnum)
                disp([syl ': EMPTY (not in any regexp motif']);
                continue
            end
            
            N_regexp=sum(~isnan(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_RegExpr.AllDays_RegExpr.baseline_MUSC.data_WithOutlier{motifnum}.FFvals(:, posinmotif))); % all samples for this syl that are not nan
                N_all=numel(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.EpochData_MUSC.Baseline.(syl).rawFF_WithinTimeWindow);
            
            if N_regexp<N_all;
                disp(['[motif' num2str(motifnum) '|pos' num2str(posinmotif) '] ' syl ': ' num2str(N_regexp) '/' num2str(N_all)]);
            else
                disp(['[motif' num2str(motifnum) '|pos' num2str(posinmotif) '] -- GOOD -- ' syl ':  ' num2str(N_regexp) '/' num2str(N_all)]);
            end
            
            
            
        end
    end
end
        
%% SAMPLE SIZES FOR ALL SYLLABLES ACROSS ALL DAYS
count=1;
SubplotsPerFig=6;
subplotrows=2;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];




for i=1:NumBirds;
    
    birdname= SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts = length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts
        
        exptname = SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        SylsUnique=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.SylFields_Unique;
        targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
        NumDays=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targsyl).FFvals);
        
            [fignums_alreadyused, hfigs, count]=lt_plot_MultSubplotsFigs(SubplotsPerFig, subplotrows, subplotcols, fignums_alreadyused, hfigs, count);
        title([birdname '-' exptname]);
        
        plotcols=lt_make_plot_colors(length(SylsUnique),0,0);
        
        hplot=[];
        for j=1:length(SylsUnique);
            
            syl=SylsUnique{j};
            
            SampleSizeAllDays=[];
            for k=1:NumDays;
                if length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals)<k
                    continue
                end
                    
                if isempty(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{k});
                    sampsize=nan;
                else
                    sampsize=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(syl).FFvals{k});
                end
                
                SampleSizeAllDays=[SampleSizeAllDays sampsize];
                
            end
            
            % === PLOT
            hplot(j)=lt_plot(1:length(SampleSizeAllDays), SampleSizeAllDays, {'LineStyle','-', 'Color', plotcols{j}});
            
        end
        
        % ===== annotate plot with lines
        Ylim=ylim;
        Yrange=Ylim(2)-Ylim(1);
        Xlim=xlim;
        
        % WN on and off
        WNonday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
        line([WNonday-0.5 WNonday-0.5], ylim, 'Color' ,'r');
        WNoffday=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOffInd;
        line([WNoffday+0.5 WNoffday+0.5],ylim,'Color','r');
        lt_plot_text(WNonday, Ylim(1), 'WNon/off','r')
        
        % day I called consol start
        try
            day1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.ConsolStartInd;
            day2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.ConsolEndInd;
            
            line([day1-0.5 day1-0.5], ylim);
            line([day2+0.5 day2+0.5], ylim);
            
            lt_plot_text(day1, Ylim(1), 'consol','b')
        catch err
        end
        
        % multidir
        if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'MultiDir_OneDayBeforeStart');
            date1_ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_OneDayBeforeStart_Ind;
            date2_ind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_LastDay_Ind;
            
            line([date1_ind+0.4 date1_ind+0.4], ylim, 'Color', 'k');
            line([date2_ind+0.4 date2_ind+0.4], ylim, 'Color', 'k');
            lt_plot_text(date1_ind, Ylim(1), 'multidir','k')
        end
        
        % LMAN inactivation days
        if  SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.LMANinactivated==1;
            targsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            muscdays_inds=find(~isnan(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(targsyl).meanFF_DevFromBase_WithinTimeWindow));
            
            % plot
            plot(muscdays_inds, Ylim(1)+Yrange/6, '^k', 'MarkerSize', 8);
            plot(muscdays_inds, Ylim(2)-Yrange/6, 'vk', 'MarkerSize', 8);
            
            %             for k=1:length(muscdays_inds);
            % %                 line([muscdays_inds(k) muscdays_inds(k)], ylim, 'Color',[0.7 0.7 0.7]);
            %                 line([muscdays_inds(k) muscdays_inds(k)], ylim, 'LineStyle','--','Color',[0.3 0.3 0.3]);
            %             end
            
            lt_plot_text(1, Ylim(1)+Yrange/8, 'arrowhead: MUSC', 'k');
            
                        
        end
        
        % annotate target learning direction (as detected automatically)
        targ_learn_dir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
        
        lt_plot_text(2, Ylim(2)-Yrange/8, ['targ learn dir: ' num2str(targ_learn_dir)], 'r');
        
        
        
        % ========= ANNOTATE 50 syls
        line(xlim, [50 50],'Color','k');
        
        legend(hplot, SylsUnique);
    end
end
    


%% [DO NOT USE!!!] REMOVE EARLY DAYS FROM EXPERIMENTS THAT LACK SONGS, SHOW NON-LEARNING EFFECT - DOES NOT WORK!! REMOVES THOSE DAYS INSTEAD OF REPLACING WITH NAN OR EMPTY
% NOTE: ALSO REMOVES FIELDS FROM STRUCTURES THAT MIGHT BE CONTRADICTORY,
% THAT I NEVER USE.

PARAMS.EarlyDays.threshold=0.8; % zscore, in opposite direction
SeqDepPitch_AcrossBirds_orig=SeqDepPitch_AcrossBirds;
[SeqDepPitch_AcrossBirds, PARAMS] = lt_seq_dep_pitch_ACROSSBIRDS_EarlyDays(SeqDepPitch_AcrossBirds, PARAMS);



%% SEARCH FOR SYLS that pass criteria

lt_seq_dep_pitch_ACROSSBIRDS_FindSyls; % go in and modify


%% ===== distribution of time windows used to calcualte pitch


NumBirds=length(SeqDepPitch_AcrossBirds.birds);
WindSizeAll=[];

for i=1:NumBirds
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    birdname = SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    for ii=1:numexpts
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        numtimewindow=length(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.pc_time_window_list{1}); % on day one, for all seq and syl
        
        tbin=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.DayRawDat.pc_T{1}(2) ...
            - SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.DayRawDat.pc_T{1}(1);
        
        % tbin should be 0.000125s
        if (0)
            disp(tbin);
        end
        
        for j=1:numtimewindow
            
            origsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.OrigNoteID{1}(j);
            origsyl=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.DayRawDat.syllables{origsyl}; % string
            
            % skip if this syl is a noisy syl (i.e. j, n, ...)
%             if strcmp(origsyl, 'j') | strcmp(origsyl, 'n') | strcmp(origsyl, 'k') | strcmp(origsyl, 'l')
%                 continue
%             end
%             
            
            timewind=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.SeqFilter.pc_time_window_list{1}(:,j);
            
            
            timediff=timewind(2)-timewind(1)+1; % bins
            timediff=timediff*tbin*1000; % msec
            
             if timediff<10
                disp([birdname '-' exptname '-' origsyl]);
            end
           WindSizeAll=[WindSizeAll timediff];
            
            
            
        end
        
    end
end

% === plot
lt_figure; lt_plot_histogram(WindSizeAll)

lt_plot_annotation(1, ['mean=' num2str(mean(WindSizeAll)) '; std=' num2str(std(WindSizeAll))], 'r')

%% PLOT STATS ABOUT THE EXPERIMENTS - E.G. SAMPLE SIZES, DAYS WHEN THINGS HAPPENED, ETC
close all;
[SeqDepPitch_AcrossBirds,PARAMS]= lt_seq_dep_pitch_ACROSSBIRDS_METASTATS(SeqDepPitch_AcrossBirds, PARAMS);

close all;
[SeqDepPitch_AcrossBirds_filtered,PARAMS]= lt_seq_dep_pitch_ACROSSBIRDS_METASTATS(SeqDepPitch_AcrossBirds_filtered, PARAMS);




%% PLOT SPECTROGRAMS OF MOTIFS/SONGS


lt_seq_dep_pitch_ACROSSBIRDS_PlotSpectrogram % IN PROGRESS



