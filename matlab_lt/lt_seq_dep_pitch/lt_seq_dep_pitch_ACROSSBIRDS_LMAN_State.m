%% LT 2/29/16 - using output structures from the multidir and samedir analyses
function lt_seq_dep_pitch_ACROSSBIRDS_LMAN_State(SeqDepPitch_AcrossBirds, Params)


NumBirds=length(SeqDepPitch_AcrossBirds.birds);

normToSingleDirLastDay=0; % difference from today

%% Extract data for desired experiments
% ---- criteria:
% 1) single --> bidir transitions
% 2) 2 targets are same-type

figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

count=1;

% ===== COLLECT DAY TO DAY FLUCTUATION ACROSS ALL EXPERIMENTS
Targ1_Diff_MUSC=[]; % -- in targ 1, diff from day a to day b, in MP expressed pitch.
Targ2_Diff_MUSC=[];

Targ1_Diff_AFP=[];
Targ2_Diff_AFP=[];

Dur_Day_Diff=[];


for i=1:NumBirds
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    numexpts=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    for ii=1:numexpts
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        % ======= PASSES CRITERIA?
        %         bidirPreStartInd=[];
        %         samedirPreStartInd=[];
        %
        %         % -- collect dates
        %         if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'MultiDir_OneDayBeforeStart_Ind');
        %             bidirPreStartInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_OneDayBeforeStart_Ind;
        %         end
        %
        %         if isfield(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES, 'SameDir_OneDayBeforeStart_Ind');
        %             samedirPreStartInd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.SameDir_OneDayBeforeStart_Ind;
        %         end
        %
        %         % --- skip if no bidir
        %         if isempty(bidirPreStartInd)
        %             continue
        %         end
        
        % -- continue only if starts with single dir, then to bidir (no
        % intervening samedir)
        numtargs=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.NumTargs;
        if numtargs==0 || numtargs==4
            % then is good
            disp(['GOOD: ' birdname '-' exptname]);
            
            % --- make sure only two bidir syls
            bidirSyls=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.MultiDirSyls;
            if length(bidirSyls)>2
               disp(['DROPPED (TOO MANY BIDIR SYLS): ' birdname '-' exptname]); 
                continue
            end
            
            % --- collect targets
            targ1=bidirSyls{1};
            targ2=bidirSyls{2};
            % - which one is orig targ?
            targorig=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targsyl;
            if strcmp(targ1, targorig)==1
                % then keep same
            elseif strcmp(targ2, targorig)==1
               % -- switch
               tmp=targ1;
               targ1=targ2;
               targ2=tmp;
            else
                disp('SOMETHING WRONG////...');
            end

            % --- make sure they are same type
            if SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Syl_ID_Dimensions.(targ2).similar_to_targ==0
               disp(['DROPPED (diff syls) ' birdname '-' exptname]); 
                continue
            end
            
            
            
            % ============================= READY TO COLLECT
            WNday1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.WNTimeOnInd;
            BidirDay1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_OneDayBeforeStart_Ind+1;
            BidirDayEnd=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.DATES.MultiDir_LastDay_Ind;
            
            DaysWithMusc=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.Params.PlotLearning.MUSCIMOL.DaysWithInact_inds;
            
            % ======= FFVALS
            % -- check
            tmp=find(~isnan(SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(targ1).meanFF_DevFromBase_WithinTimeWindow));
            assert(sum(DaysWithMusc-tmp)==0, 'sdafasdfioasdhiof');
            
            % --- MUSC
            FF_musc_targ1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(targ1).meanFF_DevFromBase_WithinTimeWindow(DaysWithMusc);
            FF_musc_targ2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix_MUSC.(targ2).meanFF_DevFromBase_WithinTimeWindow(DaysWithMusc);
            
            % --- PBS
            FF_pbs_targ1=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targ1).meanFF_DevFromBase_WithinTimeWindow(DaysWithMusc);
            FF_pbs_targ2=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.Data_PlotLearning.AllDays_PlotLearning.DataMatrix.(targ2).meanFF_DevFromBase_WithinTimeWindow(DaysWithMusc);
            
            % --- invert if learning sign neg
            targdir=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.INFORMATION.targ_learn_dir;
            
            FF_musc_targ1=FF_musc_targ1*targdir;
            FF_musc_targ2=FF_musc_targ2*targdir;
            FF_pbs_targ1=FF_pbs_targ1*targdir;
            FF_pbs_targ2=FF_pbs_targ2*targdir;
            
            % --- get AFP
            FF_afp_targ1=FF_pbs_targ1-FF_musc_targ1;
            FF_afp_targ2=FF_pbs_targ2-FF_musc_targ2;
            

            % ================= COLLECT RAW DAT ACROSS EXPERIMENTS
            DATSTRUCT.data(count).WNday1=WNday1;
            DATSTRUCT.data(count).BidirDay1=BidirDay1;
            DATSTRUCT.data(count).BidirDayEnd=BidirDayEnd;
            DATSTRUCT.data(count).FF_musc_targ1=FF_musc_targ1;
            DATSTRUCT.data(count).FF_musc_targ2=FF_musc_targ2;
            DATSTRUCT.data(count).FF_pbs_targ1=FF_pbs_targ1;
            DATSTRUCT.data(count).FF_pbs_targ2=FF_pbs_targ2;
            DATSTRUCT.data(count).DaysWithMusc=DaysWithMusc;
            
            % ======================= PLOT FOR THIS EXPT
            
            % only take days starting from last single dir day to end of
            % bidir.
            inds=find(DaysWithMusc>=BidirDay1 & DaysWithMusc<=BidirDayEnd);
            inds=[inds(1)-1 inds]; % adding the last day of single dir;

            % --- MUSC
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title([birdname(1:4) '-' exptname(end-4:end) '-MUSC']);
            xlabel('ctxt1 (in this dir)');
            ylabel('ctxt2');            
            Xtarg1=FF_musc_targ1(inds);
            Xtarg2=FF_musc_targ2(inds);
            color='r';
            
            if normToSingleDirLastDay==1;
               Xtarg1=Xtarg1-Xtarg1(1);
               Xtarg2=Xtarg2-Xtarg2(1);
            end
            
            
            plot(Xtarg1(1), Xtarg2(1), 'o-k'); % - plot last single dir day
            plot(Xtarg1(1:2), Xtarg2(1:2), '-k'); % - plot last single dir day
            plot(Xtarg1(2:end), Xtarg2(2:end), 'o-', 'Color',color); % plot bidir days;
            lt_plot_zeroline;
            lt_plot_zeroline_vert;
            
            % collect
            Targ1_Diff_MUSC=[Targ1_Diff_MUSC diff(Xtarg1)]; % -- in targ 1, diff from day a to day b, in MP expressed pitch.
            Targ2_Diff_MUSC=[Targ2_Diff_MUSC diff(Xtarg2)];


            % --- AFP
            Xtarg1=FF_afp_targ1(inds);
            Xtarg2=FF_afp_targ2(inds);
            color='b';
                
                        if normToSingleDirLastDay==1;
               Xtarg1=Xtarg1-Xtarg1(1);
               Xtarg2=Xtarg2-Xtarg2(1);
                        end

            
            plot(Xtarg1(1), Xtarg2(1), 'o-k'); % - plot last single dir day
            plot(Xtarg1(1:2), Xtarg2(1:2), '-k'); % - plot last single dir day
            plot(Xtarg1(2:end), Xtarg2(2:end), 'o-', 'Color',color); % plot bidir days;
            lt_plot_zeroline;
            lt_plot_zeroline_vert;
            
            
            Targ1_Diff_AFP=[Targ1_Diff_AFP diff(Xtarg1)];
            Targ2_Diff_AFP=[Targ2_Diff_AFP diff(Xtarg2)];

            % ------------- COLLECT OTHER STUFF
            Dur_Day_Diff=[Dur_Day_Diff diff(DaysWithMusc(inds))]; % day diff
            
            count=count+1;
            
        end
    end
    
end

%% ==== Plot correaltion of inactivationDay to inactivationDay differences

lt_figure; hold on;
title('AFP');
xlabel('first targ, diff (hz)');
ylabel('second targ');
X=Targ1_Diff_AFP;
Y=Targ2_Diff_AFP;
lt_regress(Y, X, 1, 0, 1, 1, 'b');

lt_figure; hold on;
X=Targ1_Diff_MUSC;
Y=Targ2_Diff_MUSC;
lt_regress(Y, X, 1, 0, 1, 1, 'r');


%%  ===== PLOT FIRST AND SECOND TARG AND END OF SINGLE DIR AND START OF BIDIR

NumExpts=length(DATSTRUCT.data);
MaxNumDaysFromEdge=6; % i.e. +/- from the transition

    DATSTRUCT2.data.SingleEnd.FF_musc_targ1=[];
    DATSTRUCT2.data.SingleEnd.FF_musc_targ2=[];
    DATSTRUCT2.data.SingleEnd.FF_pbs_targ1=[];
    DATSTRUCT2.data.SingleEnd.FF_pbs_targ2=[];
    DATSTRUCT2.data.BidirStart.FF_musc_targ1=[];
    DATSTRUCT2.data.BidirStart.FF_musc_targ2=[];
    DATSTRUCT2.data.BidirStart.FF_pbs_targ1=[];
    DATSTRUCT2.data.BidirStart.FF_pbs_targ2=[];

for i=1:NumExpts;
    
    bidirDay1=DATSTRUCT.data(i).BidirDay1;
    DaysWithMusc=DATSTRUCT.data(i).DaysWithMusc;

    % ==== END of single dir
    window=[bidirDay1-MaxNumDaysFromEdge bidirDay1-1];
    inds=find(DaysWithMusc>=window(1) & DaysWithMusc<=window(2));
    
    % - extract data
    ffmusc1=DATSTRUCT.data(i).FF_musc_targ1(inds);
    ffmusc2=DATSTRUCT.data(i).FF_musc_targ2(inds);
    ffpbs1=DATSTRUCT.data(i).FF_pbs_targ1(inds);
    ffpbs2=DATSTRUCT.data(i).FF_pbs_targ2(inds);
    
    DATSTRUCT2.data.SingleEnd.FF_musc_targ1=[DATSTRUCT2.data.SingleEnd.FF_musc_targ1 mean(ffmusc1)];
    DATSTRUCT2.data.SingleEnd.FF_musc_targ2=[DATSTRUCT2.data.SingleEnd.FF_musc_targ2 mean(ffmusc2)];
    DATSTRUCT2.data.SingleEnd.FF_pbs_targ1=[DATSTRUCT2.data.SingleEnd.FF_pbs_targ1 mean(ffpbs1)];
    DATSTRUCT2.data.SingleEnd.FF_pbs_targ2=[DATSTRUCT2.data.SingleEnd.FF_pbs_targ2 mean(ffpbs2)];
    
    
    % ==== START bidir
    window=[bidirDay1 bidirDay1+MaxNumDaysFromEdge-1];
    inds=find(DaysWithMusc>=window(1) & DaysWithMusc<=window(2));
    
    % - extract data
    ffmusc1=DATSTRUCT.data(i).FF_musc_targ1(inds);
    ffmusc2=DATSTRUCT.data(i).FF_musc_targ2(inds);
    ffpbs1=DATSTRUCT.data(i).FF_pbs_targ1(inds);
    ffpbs2=DATSTRUCT.data(i).FF_pbs_targ2(inds);
    
    DATSTRUCT2.data.BidirStart.FF_musc_targ1=[DATSTRUCT2.data.BidirStart.FF_musc_targ1 mean(ffmusc1)];
    DATSTRUCT2.data.BidirStart.FF_musc_targ2=[DATSTRUCT2.data.BidirStart.FF_musc_targ2 mean(ffmusc2)];
    DATSTRUCT2.data.BidirStart.FF_pbs_targ1=[DATSTRUCT2.data.BidirStart.FF_pbs_targ1 mean(ffpbs1)];
    DATSTRUCT2.data.BidirStart.FF_pbs_targ2=[DATSTRUCT2.data.BidirStart.FF_pbs_targ2 mean(ffpbs2)];
    
end

% === PLOT
lt_figure; hold on;
xlabel('end single dir -- start bidir');

% --- END SINGLE DIR
% targ 1
x=1;
ffpbs=DATSTRUCT2.data.SingleEnd.FF_pbs_targ1;
ffmusc=DATSTRUCT2.data.SingleEnd.FF_musc_targ1;

plot([x-0.1 x+0.1], [ffpbs' ffmusc'], 'ok-');
lt_plot_bar(x-0.1, mean(ffpbs), {'Errors', lt_sem(ffpbs), 'Color','k', 'BarWidth', 0.2});
lt_plot_bar(x+0.1, mean(ffmusc), {'Errors', lt_sem(ffmusc), 'Color','r','BarWidth', 0.2});

% targ 1
x=2;
ffpbs=DATSTRUCT2.data.SingleEnd.FF_pbs_targ2;
ffmusc=DATSTRUCT2.data.SingleEnd.FF_musc_targ2;

plot([x-0.1 x+0.1], [ffpbs' ffmusc'], 'ok-');
lt_plot_bar(x-0.1, mean(ffpbs), {'Errors', lt_sem(ffpbs), 'Color','k', 'BarWidth', 0.2});
lt_plot_bar(x+0.1, mean(ffmusc), {'Errors', lt_sem(ffmusc), 'Color','r','BarWidth', 0.2});

% --- START BIDIR
% targ 1
x=4;
ffpbs=DATSTRUCT2.data.BidirStart.FF_pbs_targ1;
ffmusc=DATSTRUCT2.data.BidirStart.FF_musc_targ1;

plot([x-0.1 x+0.1], [ffpbs' ffmusc'], 'ok-');
lt_plot_bar(x-0.1, mean(ffpbs), {'Errors', lt_sem(ffpbs), 'Color','k', 'BarWidth', 0.2});
lt_plot_bar(x+0.1, mean(ffmusc), {'Errors', lt_sem(ffmusc), 'Color','r','BarWidth', 0.2});

% targ 1
x=5;
ffpbs=DATSTRUCT2.data.BidirStart.FF_pbs_targ2;
ffmusc=DATSTRUCT2.data.BidirStart.FF_musc_targ2;

plot([x-0.1 x+0.1], [ffpbs' ffmusc'], 'ok-');
lt_plot_bar(x-0.1, mean(ffpbs), {'Errors', lt_sem(ffpbs), 'Color','k', 'BarWidth', 0.2});
lt_plot_bar(x+0.1, mean(ffmusc), {'Errors', lt_sem(ffmusc), 'Color','r','BarWidth', 0.2});


% ---- diff across phases?
ffmusc1=DATSTRUCT2.data.SingleEnd.FF_musc_targ1; % targ 1, musc (single)
ffmusc2=DATSTRUCT2.data.BidirStart.FF_musc_targ1; % (bidir)

p=signrank(ffmusc1, ffmusc2);
line([1.1 4.1], [max(ffmusc1)+50 max(ffmusc1)+50], 'Color','r');
lt_plot_text(2.3, [max(ffmusc1)+50]+10, ['p=' num2str(p)],'r');


% ---- diff across phases?
ffmusc1=DATSTRUCT2.data.SingleEnd.FF_musc_targ2; % targ 1, musc (single)
ffmusc2=DATSTRUCT2.data.BidirStart.FF_musc_targ2; % (bidir)

p=signrank(ffmusc1, ffmusc2);
line([2.1 5.1], [max(ffmusc1)+50 max(ffmusc1)+50], 'Color','m');
lt_plot_text(2.3, [max(ffmusc1)+50]+10, ['p=' num2str(p)],'m');

lt_plot_annotation(1, ['using ' num2str(MaxNumDaysFromEdge) ' days +/- from single --> bidir switch [mean of days]'], 'b');
