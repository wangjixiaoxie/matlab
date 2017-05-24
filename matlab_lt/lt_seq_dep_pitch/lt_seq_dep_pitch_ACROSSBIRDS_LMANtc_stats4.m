%% LT 2/29/16 - using output structures from the multidir and samedir analyses
function lt_seq_dep_pitch_ACROSSBIRDS_LMANtc_stats4(SeqDepPitch_AcrossBirds,OUTPUT_multidir, OUTPUT_samedir, OUTPUT_singleTarg, DATSTRUCT_multidir, DATSTRUCT_samedir, DATSTRUCT_singleTarg,  DayBinToUse, AdHocExpt, TakeAvg, ExpToAvg)

figcount=1;
subplotrows=4;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

% DayBinToUse=3;

if ~exist('AdHocExpt', 'var');
    AdHocExpt={};
end

if ~exist('TakeAvg', 'var');
    TakeAvg={};
end


%%
lt_figure;

lt_plot_annotation(1, 'note: this contains all experiments, even if multiple of same time for same Bird and Expt (hand entered)', 'k');
lt_plot_annotation(2, 'wh25 last bidir expt is missing hit rate!!', 'k');
lt_plot_annotation(3, 'paired lines only go between adjacent datapoints - make them skip!!', 'k');



%% ====== GET TIMECOURSE FOR EACH EXPERIMENT (i.e. combine single targ, bidir, and samedir


% FIRST GET LIST OF ALL POTENTIAL EXPERIMENTS ACROSS ALL THINGS
NumTotalExpts=0;
% Birdnames={};
% Exptnames={};
GLOB_STRUCT=struct; % INITIATE OUTPUT STRUCTURE

for i=1:length(SeqDepPitch_AcrossBirds.birds);
    numexpt=length(SeqDepPitch_AcrossBirds.birds{i}.experiment);
    
    birdname=SeqDepPitch_AcrossBirds.birds{i}.birdname;
    
    for ii=1:numexpt
        exptname=SeqDepPitch_AcrossBirds.birds{i}.experiment{ii}.ExptID;
        
        % --- COLLECT
        NumTotalExpts=NumTotalExpts+1;
        %         Birdnames=[Birdnames birdname];
        %         Exptnames=[Exptnames exptname];
        
        GLOB_STRUCT.info(NumTotalExpts).birdname=birdname;
        GLOB_STRUCT.info(NumTotalExpts).exptname=exptname;
        
        
    end
end

% initiate so have maximum length
GLOB_STRUCT.data.multiDir(NumTotalExpts).consolidation=[];
GLOB_STRUCT.data.singleTarg(NumTotalExpts).consolidation=[];
GLOB_STRUCT.data.sameDir(NumTotalExpts).consolidation=[];



% ++++++++++ SECOND, go through each of the outputs structures, and slot it
% into the correct bird/expt

% ============== 1) single target
exptnums=OUTPUT_singleTarg.INFORMATION(DayBinToUse).experimentNum;

for i=1:length(exptnums);
    
    exptnum=exptnums(i);
    
    % -- which bird expt?
    birdname=DATSTRUCT_singleTarg.INFORMATION(exptnum).birdname;
    exptname=DATSTRUCT_singleTarg.INFORMATION(exptnum).exptname;
    
    % -- what day is this period relative to entire experiment
    day1_relWN=DATSTRUCT_singleTarg.INFORMATION(exptnum).consolPeriod(1);
    
    % -- collect data
    ffPBS=OUTPUT_singleTarg.firsttarget.FFrelBaseInBin(DayBinToUse).PBS(i);
    ffMUSC=OUTPUT_singleTarg.firsttarget.FFrelBaseInBin(DayBinToUse).MUSC(i);
    
    consolidation=ffMUSC/ffPBS;
    
    hitrate_PBS=OUTPUT_singleTarg.firsttarget.otherstats(DayBinToUse).hitratePBS(i);
    hitrate_MUSC=OUTPUT_singleTarg.firsttarget.otherstats(DayBinToUse).hitrateMUSC(i);
    
    % ---- SAVE IT TO CORRECT SLOT
    
    for j=1:length(GLOB_STRUCT.info)
        if strcmp(GLOB_STRUCT.info(j).birdname, birdname) ...
                & strcmp(GLOB_STRUCT.info(j).exptname, exptname);
            
            % -- then this is the correct place to put it
            GLOB_STRUCT.data.singleTarg(j).ffPBS=ffPBS;
            GLOB_STRUCT.data.singleTarg(j).ffMUSC=ffMUSC;
            GLOB_STRUCT.data.singleTarg(j).consolidation=consolidation;
            GLOB_STRUCT.data.singleTarg(j).consolidation_day1_relWN=day1_relWN;
            GLOB_STRUCT.data.singleTarg(j).exptnum_WithinSingleTarg=exptnum;
            GLOB_STRUCT.data.singleTarg(j).birdname=birdname;
            GLOB_STRUCT.data.singleTarg(j).exptname=exptname;
            
            GLOB_STRUCT.data.singleTarg(j).hrate_PBS=hitrate_PBS;
            GLOB_STRUCT.data.singleTarg(j).hrate_MUSC=hitrate_MUSC;
            
        end
    end
end


% ============== 1) twotarg, bidir
exptnums=OUTPUT_multidir.INFORMATION(DayBinToUse).experimentNum;

for i=1:length(exptnums);
    
    exptnum=exptnums(i);
    
    % -- which bird expt?
    birdname=DATSTRUCT_multidir.INFORMATION(exptnum).birdname;
    exptname=DATSTRUCT_multidir.INFORMATION(exptnum).exptname;
    
    % -- what day is this period relative to entire experiment
    phase_day1_relWN=DATSTRUCT_multidir.INFORMATION(exptnum).day1_FromStartWN;
    
    % -- collect data
    ff_combined_PBS=OUTPUT_multidir.separation(DayBinToUse).PBS(i);
    ff_combined_MUSC=OUTPUT_multidir.separation(DayBinToUse).MUSC(i);
    
    ff_first_PBS=OUTPUT_multidir.firsttarget.FFrelBaseInBin(DayBinToUse).PBS(i);
    ff_first_MUSC=OUTPUT_multidir.firsttarget.FFrelBaseInBin(DayBinToUse).MUSC(i);
    ff_second_PBS=OUTPUT_multidir.secondtarg.FFrelBaseInBin(DayBinToUse).PBS(i);
    ff_second_MUSC=OUTPUT_multidir.secondtarg.FFrelBaseInBin(DayBinToUse).MUSC(i);
    
    consol_first=ff_first_MUSC/ff_first_PBS;
    consol_second=ff_second_MUSC/ff_second_PBS;
    
    consolidation=ff_combined_MUSC/ff_combined_PBS;
    
    hitrate_comb_PBS=OUTPUT_multidir.combined.otherstats(DayBinToUse).hitratePBS(i);
    hitrate_comb_MUSC=OUTPUT_multidir.combined.otherstats(DayBinToUse).hitrateMUSC(i);
    
    % ---- SAVE IT TO CORRECT SLOT
    
    for j=1:length(GLOB_STRUCT.info)
        if strcmp(GLOB_STRUCT.info(j).birdname, birdname) ...
                & strcmp(GLOB_STRUCT.info(j).exptname, exptname);
            
            % -- then this is the correct place to put it
            GLOB_STRUCT.data.multiDir(j).ff_combined_PBS=ff_combined_PBS;
            GLOB_STRUCT.data.multiDir(j).ff_combined_MUSC=ff_combined_MUSC;
            GLOB_STRUCT.data.multiDir(j).consolidation=consolidation;
            GLOB_STRUCT.data.multiDir(j).phase_day1_relWN=phase_day1_relWN;
            GLOB_STRUCT.data.multiDir(j).exptnum_WithinSingleTarg=exptnum;
            GLOB_STRUCT.data.multiDir(j).birdname=birdname;
            GLOB_STRUCT.data.multiDir(j).exptname=exptname;
            
            GLOB_STRUCT.data.multiDir(j).hrate_combined_PBS=hitrate_comb_PBS;
            GLOB_STRUCT.data.multiDir(j).hrate_combined_MUSC=hitrate_comb_MUSC;
            
            GLOB_STRUCT.data.multiDir(j).ff_first_PBS=ff_first_PBS;
            GLOB_STRUCT.data.multiDir(j).ff_first_MUSC=ff_first_MUSC;
            GLOB_STRUCT.data.multiDir(j).ff_second_PBS=ff_second_PBS;
            GLOB_STRUCT.data.multiDir(j).ff_second_MUSC=ff_second_MUSC;
            GLOB_STRUCT.data.multiDir(j).consol_first=consol_first;
            GLOB_STRUCT.data.multiDir(j).consol_second=consol_second;
            
        end
    end
    
end


% ============== 1) twotarg, samedir
exptnums=OUTPUT_samedir.INFORMATION(DayBinToUse).experimentNum;

for i=1:length(exptnums);
    
    exptnum=exptnums(i);
    
    % -- which bird expt?
    birdname=DATSTRUCT_samedir.INFORMATION(exptnum).birdname;
    exptname=DATSTRUCT_samedir.INFORMATION(exptnum).exptname;
    
    % -- what day is this period relative to entire experiment
    phase_day1_relWN=DATSTRUCT_samedir.INFORMATION(exptnum).day1_FromStartWN;
    
    % -- collect data
    ff_combined_PBS=OUTPUT_samedir.separation(DayBinToUse).PBS(i);
    ff_combined_MUSC=OUTPUT_samedir.separation(DayBinToUse).MUSC(i);
    
    ff_first_PBS=OUTPUT_samedir.firsttarget.FFrelBaseInBin(DayBinToUse).PBS(i);
    ff_first_MUSC=OUTPUT_samedir.firsttarget.FFrelBaseInBin(DayBinToUse).MUSC(i);
    ff_second_PBS=OUTPUT_samedir.secondtarg.FFrelBaseInBin(DayBinToUse).PBS(i);
    ff_second_MUSC=OUTPUT_samedir.secondtarg.FFrelBaseInBin(DayBinToUse).MUSC(i);
    
    consol_first=ff_first_MUSC/ff_first_PBS;
    consol_second=ff_second_MUSC/ff_second_PBS;
    
    consolidation=ff_combined_MUSC/ff_combined_PBS;
    
    hitrate_comb_PBS=OUTPUT_samedir.combined.otherstats(DayBinToUse).hitratePBS(i);
    hitrate_comb_MUSC=OUTPUT_samedir.combined.otherstats(DayBinToUse).hitrateMUSC(i);
    
    % ---- SAVE IT TO CORRECT SLOT
    
    for j=1:length(GLOB_STRUCT.info)
        if strcmp(GLOB_STRUCT.info(j).birdname, birdname) ...
                & strcmp(GLOB_STRUCT.info(j).exptname, exptname);
            
            % -- then this is the correct place to put it
            GLOB_STRUCT.data.sameDir(j).ff_combined_PBS=ff_combined_PBS;
            GLOB_STRUCT.data.sameDir(j).ff_combined_MUSC=ff_combined_MUSC;
            GLOB_STRUCT.data.sameDir(j).consolidation=consolidation;
            GLOB_STRUCT.data.sameDir(j).phase_day1_relWN=phase_day1_relWN;
            GLOB_STRUCT.data.sameDir(j).exptnum_WithinSingleTarg=exptnum;
            
            GLOB_STRUCT.data.sameDir(j).hrate_combined_PBS=hitrate_comb_PBS;
            GLOB_STRUCT.data.sameDir(j).hrate_combined_MUSC=hitrate_comb_MUSC;
            
            GLOB_STRUCT.data.sameDir(j).birdname=birdname;
            GLOB_STRUCT.data.sameDir(j).exptname=exptname;
            
            GLOB_STRUCT.data.sameDir(j).ff_first_PBS=ff_first_PBS;
            GLOB_STRUCT.data.sameDir(j).ff_first_MUSC=ff_first_MUSC;
            GLOB_STRUCT.data.sameDir(j).ff_second_PBS=ff_second_PBS;
            GLOB_STRUCT.data.sameDir(j).ff_second_MUSC=ff_second_MUSC;
            GLOB_STRUCT.data.sameDir(j).consol_first=consol_first;
            GLOB_STRUCT.data.sameDir(j).consol_second=consol_second;
            
        end
    end
    
    
end

% ======== 2) ADD THINGS BY HAND - i.e. does not fit into previous analyses
% (e.g. two phases of the same thing in one experiment

numtoadd=length(AdHocExpt)/5;
if mod(length(AdHocExpt), 5)~=0
    disp('PROBLEM!');
    asdcae4er;
end

for i=1:numtoadd;
    % extract info
    birdname=AdHocExpt{5*i-4};
    exptname=AdHocExpt{5*i-3};
    expttype=AdHocExpt{5*i-2};
    varsToAdd=AdHocExpt{5*i-1};
    NoteToSelf=AdHocExpt{5*i};
    
    % --- find correct slot
    for j=1:length(GLOB_STRUCT.info)
        
        % get birdname and exptname for this expt
        if strcmp(GLOB_STRUCT.info(j).birdname, birdname) & ...
                strcmp(GLOB_STRUCT.info(j).exptname, exptname)
            
            % === YES! - insert into this structure
            
            for k=1:length(varsToAdd)/2;
                
                GLOB_STRUCT.data.(expttype)(j).(varsToAdd{2*k-1})=varsToAdd{2*k};
            end
            
            % calculate consolidation
            try GLOB_STRUCT.data.(expttype)(j).consolidation
            catch err
                GLOB_STRUCT.data.(expttype)(j).consolidation= ...
                    GLOB_STRUCT.data.(expttype)(j).ff_combined_MUSC/...
                    GLOB_STRUCT.data.(expttype)(j).ff_combined_PBS;
            end
            
            GLOB_STRUCT.data.(expttype)(j).birdname=birdname;
            GLOB_STRUCT.data.(expttype)(j).exptname=exptname;
            
            % ---
            GLOB_STRUCT.data.(expttype)(j).consol_first=...
                GLOB_STRUCT.data.(expttype)(j).ff_first_MUSC/...
                GLOB_STRUCT.data.(expttype)(j).ff_first_PBS;
            
            GLOB_STRUCT.data.(expttype)(j).consol_second=...
                GLOB_STRUCT.data.(expttype)(j).ff_second_MUSC/...
                GLOB_STRUCT.data.(expttype)(j).ff_second_PBS;
            
            
            disp(['Added Ad hoc for ' birdname '-' exptname '; NOTE: ' NoteToSelf ]);
            
        end
    end
end


%% ====== TAKE AVERAGE OF EXPERIMENTS IF DESIRED
% NOTE: will delete the fields that are present for both experiments
% (leaving those that are not shared), and will replace those fields (the
% sahred one) with new values in the first of the old expeirments.

NumExpts=length(GLOB_STRUCT.info);

if TakeAvg==1
    NumToDo=length(ExpToAvg);
    
    for i=1:NumToDo
        TMPSTRUCT=struct;
        birdname1=ExpToAvg{i}{1};
        exptname1=ExpToAvg{i}{2};
        phase1=ExpToAvg{i}{3};
        birdname2=ExpToAvg{i}{4};
        exptname2=ExpToAvg{i}{5};
        phase2=ExpToAvg{i}{6};
        
        % ==== 1) find the first eexeriment
        exptnumber1=[];
        exptnumber2=[];
        for j=1:NumExpts
            if strcmp(GLOB_STRUCT.info(j).birdname, birdname1) & strcmp(GLOB_STRUCT.info(j).exptname, exptname1)
                % then is expt 1
                if ~isempty(exptnumber1)
                    disp('WHAT!!');
                end % troubleshooting
                exptnumber1=j;
            end
            
            if strcmp(GLOB_STRUCT.info(j).birdname, birdname2) & strcmp(GLOB_STRUCT.info(j).exptname, exptname2)
                if ~isempty(exptnumber2)
                    disp('WHAT!!');
                end % troubleshooting
                exptnumber2=j;
            end
            
        end
        
        % ======== 2) TAKE AVERAGE OF EXPERIMENTS
        % --- get list of things to average
        datafields=fieldnames(GLOB_STRUCT.data.(phase1));
        datafields_tmp={};
        % -- only keep datafields that are arrays
        for k=1:length(datafields);
            if isnumeric(GLOB_STRUCT.data.(phase1)(exptnumber1).(datafields{k}))
                % then keep
                datafields_tmp=[datafields_tmp datafields{k}];
            end
        end
        datafields=datafields_tmp;
        
        % --- TAKE AVERAGE BETWEEN TWO EXPERIMENTS
        for k=1:length(datafields)
            try % as sometimes one expt is missing a datafield
                TMPSTRUCT.(datafields{k})= mean([GLOB_STRUCT.data.(phase1)(exptnumber1).(datafields{k}) ...
            GLOB_STRUCT.data.(phase2)(exptnumber2).(datafields{k})]);
        
%         % --- erase the current values
%         GLOB_STRUCT.data.(phase1)(exptnumber1).(datafields{k})=[];
%         GLOB_STRUCT.data.(phase2)(exptnumber2).(datafields{k})=[];
            catch err 
            end
        end
        
        % --- erase current values (saves birdname and exptname for expt1)
        for k=1:length(datafields)
        % --- erase the current values
        GLOB_STRUCT.data.(phase1)(exptnumber1).(datafields{k})=[];
        end
        datfields2_tmp=fieldnames(GLOB_STRUCT.data.(phase2)(exptnumber2));
        for k=1:length(datfields2_tmp)
        % --- erase the current values
        GLOB_STRUCT.data.(phase2)(exptnumber2).(datfields2_tmp{k})=[];
        end
        
        
        % --- recalculate ratios (as the mean of ratios would not be
        % accurate)
        TMPSTRUCT.consolidation=TMPSTRUCT.ff_combined_MUSC./TMPSTRUCT.ff_combined_PBS;
        TMPSTRUCT.consol_first=TMPSTRUCT.ff_first_MUSC./TMPSTRUCT.ff_first_PBS;
        TMPSTRUCT.consol_second=TMPSTRUCT.ff_second_MUSC./TMPSTRUCT.ff_second_PBS;

        % ======= 3) save average
        for k=1:length(datafields)
            try % as sometimes one expt is missing a datafield
                 GLOB_STRUCT.data.(phase1)(exptnumber1).(datafields{k})= ...
                     TMPSTRUCT.(datafields{k});
            catch err 
            end
        end
    end
end


%% ====== PLOT EACH EXPERIMENT [FF]
figcount=1;
subplotrows=2;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

for i=1:NumTotalExpts;
    
    birdname=GLOB_STRUCT.info(i).birdname;
    exptname=GLOB_STRUCT.info(i).exptname;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['FF; ' birdname '-' exptname]);
    xlim([0 50]);
    ylim([-0.3 1.2]);
    
    % ==== COLLECT data
    Xall=[];
    Yall=[];
    stardatesAll=[];
    
    % ======= singleTarg?
    color='k';
    if ~isempty(GLOB_STRUCT.data.singleTarg(i).consolidation)
        
        consol=GLOB_STRUCT.data.singleTarg(i).consolidation;
        startday=GLOB_STRUCT.data.singleTarg(i).consolidation_day1_relWN;
        
        lt_plot(startday, consol, {'Color',color});
        
    end
    
    % ======= bidir?
    color='r';
    if ~isempty(GLOB_STRUCT.data.multiDir(i).consolidation)
        
        consol=GLOB_STRUCT.data.multiDir(i).consolidation;
        startday=GLOB_STRUCT.data.multiDir(i).phase_day1_relWN;
        
        lt_plot(startday, consol, {'Color',color});
        
    end
    
    % ======= samedir?
    color='b';
    if ~isempty(GLOB_STRUCT.data.sameDir(i).consolidation)
        
        consol=GLOB_STRUCT.data.sameDir(i).consolidation;
        startday=GLOB_STRUCT.data.sameDir(i).phase_day1_relWN;
        
        lt_plot(startday, consol, {'Color',color});
        
    end
    
    
    % ====== second bidir phase?
    color='r';
    if isfield(GLOB_STRUCT.data, 'multiDir2')
        if ~isempty(GLOB_STRUCT.data.multiDir2(i).consolidation)
            
            consol=GLOB_STRUCT.data.multiDir2(i).consolidation;
            startday=GLOB_STRUCT.data.multiDir2(i).phase_day1_relWN;
            
            lt_plot(startday, consol, {'Color',color});
            
        end
    end
    
end


%% ====== PLOT EACH EXPERIMENT [first FF]
figcount=1;
subplotrows=2;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

for i=1:NumTotalExpts;
    
    birdname=GLOB_STRUCT.info(i).birdname;
    exptname=GLOB_STRUCT.info(i).exptname;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['1st targ consol; ' birdname '-' exptname]);
    xlim([0 50]);
    ylim([-0.3 1.2]);
    
    % ==== COLLECT data
    Xall=[];
    Yall=[];
    stardatesAll=[];
    
    % ======= singleTarg?
    color='k';
    if ~isempty(GLOB_STRUCT.data.singleTarg(i).consolidation)
        
        consol=GLOB_STRUCT.data.singleTarg(i).ffMUSC/GLOB_STRUCT.data.singleTarg(i).ffPBS;
        startday=GLOB_STRUCT.data.singleTarg(i).consolidation_day1_relWN;
        
        lt_plot(startday, consol, {'Color',color});
        
    end
    
    % ======= bidir?
    color='r';
    if ~isempty(GLOB_STRUCT.data.multiDir(i).consolidation)
        
        consol=GLOB_STRUCT.data.multiDir(i).ff_first_MUSC/GLOB_STRUCT.data.multiDir(i).ff_first_PBS;
        startday=GLOB_STRUCT.data.multiDir(i).phase_day1_relWN;
        
        lt_plot(startday, consol, {'Color',color});
        
    end
    
    % ======= samedir?
    color='b';
    if ~isempty(GLOB_STRUCT.data.sameDir(i).consolidation)
        
        consol=GLOB_STRUCT.data.sameDir(i).ff_first_MUSC/GLOB_STRUCT.data.sameDir(i).ff_first_PBS;
        startday=GLOB_STRUCT.data.sameDir(i).phase_day1_relWN;
        
        lt_plot(startday, consol, {'Color',color});
        
    end
    
    
    % ====== second bidir phase?
    color='r';
    if isfield(GLOB_STRUCT.data, 'multiDir2')
        if ~isempty(GLOB_STRUCT.data.multiDir2(i).consolidation)
            
            consol=GLOB_STRUCT.data.multiDir2(i).ff_first_MUSC/GLOB_STRUCT.data.multiDir2(i).ff_first_PBS;
            startday=GLOB_STRUCT.data.multiDir2(i).phase_day1_relWN;
            
            lt_plot(startday, consol, {'Color',color});
            
        end
    end
    
end

%% ====== PLOT EACH EXPERIMENT [HIT RATE]
figcount=1;
subplotrows=2;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];

for i=1:NumTotalExpts;
    
    birdname=GLOB_STRUCT.info(i).birdname;
    exptname=GLOB_STRUCT.info(i).exptname;
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['hit rate; ' birdname '-' exptname]);
    xlim([0 50]);
    ylim([-0.3 1.2]);
    
    
    % +++++++++++++++++++++++++++++++++ PBS
    % ======= singleTarg?
    color='k';
    if ~isempty(GLOB_STRUCT.data.singleTarg(i).hrate_PBS)
        
        hrate=GLOB_STRUCT.data.singleTarg(i).hrate_PBS;
        startday=GLOB_STRUCT.data.singleTarg(i).consolidation_day1_relWN;
        
        lt_plot(startday, hrate, {'Color',color});
        
    end
    
    % ======= bidir?
    color='r';
    if ~isempty(GLOB_STRUCT.data.multiDir(i).hrate_combined_PBS)
        
        hrate=GLOB_STRUCT.data.multiDir(i).hrate_combined_PBS;
        startday=GLOB_STRUCT.data.multiDir(i).phase_day1_relWN;
        
        lt_plot(startday, hrate, {'Color',color});
        
    end
    
    % ======= samedir?
    color='b';
    if ~isempty(GLOB_STRUCT.data.sameDir(i).hrate_combined_PBS)
        
        hrate=GLOB_STRUCT.data.sameDir(i).hrate_combined_PBS;
        startday=GLOB_STRUCT.data.sameDir(i).phase_day1_relWN;
        
        lt_plot(startday, hrate, {'Color',color});
        
    end
    
    
    % ++++++++++++++++++++++++++++ MUSC
    % ======= singleTarg?
    color='k';
    if ~isempty(GLOB_STRUCT.data.singleTarg(i).hrate_MUSC)
        
        hrate=GLOB_STRUCT.data.singleTarg(i).hrate_MUSC;
        startday=GLOB_STRUCT.data.singleTarg(i).consolidation_day1_relWN;
        
        plot(startday, hrate, 's', 'Color',color);
        
    end
    
    % ======= bidir?
    color='r';
    if ~isempty(GLOB_STRUCT.data.multiDir(i).hrate_combined_MUSC)
        
        hrate=GLOB_STRUCT.data.multiDir(i).hrate_combined_MUSC;
        startday=GLOB_STRUCT.data.multiDir(i).phase_day1_relWN;
        
        plot(startday, hrate, 's', 'Color',color);
        
    end
    
    % ======= samedir?
    color='b';
    if ~isempty(GLOB_STRUCT.data.sameDir(i).hrate_combined_MUSC)
        
        hrate=GLOB_STRUCT.data.sameDir(i).hrate_combined_MUSC;
        startday=GLOB_STRUCT.data.sameDir(i).phase_day1_relWN;
        
        plot(startday, hrate, 's', 'Color',color);
        
    end
    
    
    %     % ====== second bidir phase?
    %     color='r';
    %     if ~isempty(GLOB_STRUCT.data.multiDir2(i))
    %
    %         hrate=GLOB_STRUCT.data.multiDir2(i).consolidation;
    %         startday=GLOB_STRUCT.data.multiDir2(i).phase_day1_relWN;
    %
    %         lt_plot(startday, hrate, {'Color',color});
    %
    %     end
    
    
end

%% ===== PLOT ALL ONTO ONE PLOT [BOTH SYLS]
if (0)
    lt_figure; hold on;
    plot_text=0;
    
    % ++++ INITIATE SUBPLOTS
    
    % ==== plot 1: singe targ --> diff dir --> same dir ---> diff dir
    lt_subplot(2,1,1); hold on;
    title('single --> incongr --> congr --> incongr');
    xlim([0 5]); ylim([-0.2 1.2])
    % ==== plot 2
    lt_subplot(2,1,2); hold on;
    title('single --> congr --> incongr -->');
    xlim([0 5]); ylim([-0.2 1.2])
    
    
    SummaryCell_plot1=cell(1, 4);
    SummaryCell_plot2=cell(1, 4);
    
    
    for i=1:NumTotalExpts
        
        birdname=GLOB_STRUCT.info(i).birdname;
        exptname=GLOB_STRUCT.info(i).exptname;
        
        %     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        %     title([birdname '-' exptname]);
        %     xlim([0 50]);
        %     ylim([-0.3 1.2]);
        
        % ==== COLLECT data
        TypeAll=[]; % 0 = single; 1 = bidir; 2 = samedir
        ConsolAll=[];
        stardatesAll=[];
        
        % ======= singleTarg?
        type=0;
        if ~isempty(GLOB_STRUCT.data.singleTarg(i).consolidation)
            
            consol=GLOB_STRUCT.data.singleTarg(i).consolidation;
            startday=GLOB_STRUCT.data.singleTarg(i).consolidation_day1_relWN;
            
            % -- collect
            TypeAll=[TypeAll type];
            ConsolAll=[ConsolAll consol];
            stardatesAll=[stardatesAll startday];
            
        end
        
        % ======= bidir?
        type = 1;
        if ~isempty(GLOB_STRUCT.data.multiDir(i).consolidation)
            
            consol=GLOB_STRUCT.data.multiDir(i).consolidation;
            startday=GLOB_STRUCT.data.multiDir(i).phase_day1_relWN;
            
            % -- collect
            TypeAll=[TypeAll type];
            ConsolAll=[ConsolAll consol];
            stardatesAll=[stardatesAll startday];
        end
        
        % ======= bidir (second phase)?
        type = 1;
        if isfield(GLOB_STRUCT.data, 'multiDir2')
            
            if ~isempty(GLOB_STRUCT.data.multiDir2(i).consolidation)
                
                consol=GLOB_STRUCT.data.multiDir2(i).consolidation;
                startday=GLOB_STRUCT.data.multiDir2(i).phase_day1_relWN;
                
                % -- collect
                TypeAll=[TypeAll type];
                ConsolAll=[ConsolAll consol];
                stardatesAll=[stardatesAll startday];
            end
        end
        
        % ======= samedir?
        type = 2;
        if ~isempty(GLOB_STRUCT.data.sameDir(i).consolidation)
            
            consol=GLOB_STRUCT.data.sameDir(i).consolidation;
            startday=GLOB_STRUCT.data.sameDir(i).phase_day1_relWN;
            
            % -- collect
            TypeAll=[TypeAll type];
            ConsolAll=[ConsolAll consol];
            stardatesAll=[stardatesAll startday];
            
        end
        
        % ++++ sort by start date of epochs
        [~, inds]=sort(stardatesAll);
        
        TypeAll=TypeAll(inds);
        ConsolAll=ConsolAll(inds);
        stardatesAll=stardatesAll(inds);
        
        %     if strcmp(birdname, 'bk34bk68') & strcmp(exptname, 'SeqDepPitchLMAN3');
        %         keyboard
        %     end
        
        % ++++ determine which figure to plot in based on order of types
        plotsMade=0;
        
        % --- figure 1?
        tmp=num2str(TypeAll);
        typeallStr=regexprep(tmp,'[^\w'']','');
        
        ind=regexp('0121', typeallStr);
        if length(typeallStr)==1,
            % so does not match multiple phases
            ind=ind(1);
        end
        
        if ~isempty(ind) & ind<3;
            plotX=ind:ind+length(typeallStr)-1;
            
            subplot(2,1,1);
            
            plot(plotX, ConsolAll, 'o-');
            plotsMade=plotsMade+1;
            if plot_text==1;
                lt_plot_text(max(plotX), ConsolAll(end), [birdname(1:4) '-' exptname]);
            end
            
            % ---- COLLECT FOR MEANS
            for k=1:length(plotX);
                SummaryCell_plot1{plotX(k)}=[SummaryCell_plot1{plotX(k)} ConsolAll(k)];
            end
            
        end
        
        % --- figure 2?
        tmp=num2str(TypeAll);
        typeallStr=regexprep(tmp,'[^\w'']','');
        
        ind=regexp('0212', typeallStr);
        if length(typeallStr)==1,
            % so does not match multiple phases
            ind=ind(1);
        end
        
        if ~isempty(ind) & ind<3;
            plotX=ind:ind+length(typeallStr)-1;
            
            subplot(2,1,2);
            
            plot(plotX, ConsolAll, 'o-');
            plotsMade=plotsMade+1;
            if plot_text==1;
                lt_plot_text(max(plotX), ConsolAll(end), [birdname(1:4) '-' exptname]);
            end
            
            % ---- COLLECT FOR MEANS
            for k=1:length(plotX);
                SummaryCell_plot2{plotX(k)}=[SummaryCell_plot2{plotX(k)} ConsolAll(k)];
            end
            
        end
        
        
        %     % ++++++++++++++++++
        %     [tmp, plotX]=ismember(TypeAll, [0 1 2 1]);
        %     if ~isempty(plotX) & issorted(plotX) & ~any(diff(plotX)-1) & plotX(1)<3
        %         %     [~, ~, plotX]=intersect(TypeAll, [0 1 2 1]);
        %         %     if ~isempty(plotX) & issorted(plotX) & ~any(diff(plotX)-1) & plotX(1)<3
        %         % means this expt follows that order. <3 is to make sure matches
        %         % the start of the cycle (so not redundant across plots). issorted
        %         % is making sure is in chronological order. any...diff.. is to make
        %         % sure they are adjacent (e.g. 1 2 3, not 1 3)
        %
        %         subplot(2,1,1);
        %
        %         plot(plotX, ConsolAll, 'o-');
        %         plotsMade=plotsMade+1;
        %         if plot_text==1;
        %             lt_plot_text(max(plotX), ConsolAll(end), [birdname(1:4) '-' exptname]);
        %         end
        %
        %         % ---- COLLECT FOR MEANS
        %         for k=1:length(plotX);
        %             SummaryCell_plot1{plotX(k)}=[SummaryCell_plot1{plotX(k)} ConsolAll(k)];
        %         end
        %
        %     end
        %
        %     % --- figure 2?
        %     [tmp, plotX]=ismember(TypeAll, [0 2 1 2]);
        %     if ~isempty(plotX) & issorted(plotX) & ~any(diff(plotX)-1) & plotX(1)<3
        %         %     [~, ~, plotX]=intersect(TypeAll, [0 2 1 2]);
        %         %     if ~isempty(plotX) & issorted(plotX) & ~any(diff(plotX)-1) & plotX(1)<3
        %         % means this expt follows that order. <3 is to make sure matches
        %         % the start of the cycle (so not redundant across plots)
        %
        %         subplot(2,1,2);
        %
        %         plot(plotX, ConsolAll, 'o-');
        %         plotsMade=plotsMade+1;
        %
        %         if plot_text==1;
        %             lt_plot_text(max(plotX), ConsolAll(end), [birdname(1:4) '-' exptname]);
        %         end
        %
        %
        %         % ---- COLLECT FOR MEANS
        %         for k=1:length(plotX);
        %             SummaryCell_plot2{plotX(k)}=[SummaryCell_plot2{plotX(k)} ConsolAll(k)];
        %         end
        %
        %     end
        % +++++++++++++++++++++++++++++++++++++++++++
        
        if plotsMade>1 & TypeAll~=0
            disp('PROBLEM - plotted >1 time, and was not simply single-targ plots');
        end
        
        if plotsMade<1
            disp(['PROBLEM!! - ' birdname '-' exptname ' failed to find correct subplot']);
        end
        
    end
    
    
    % ===== OVERLAY BARS FOR MEANS
    subplot(2,1,1);
    for i=1:length(SummaryCell_plot1);
        yvals=SummaryCell_plot1{i};
        ymean=mean(yvals);
        ysem=lt_sem(yvals);
        
        lt_plot_bar(i, ymean, {'Errors', ysem});
    end
    
    subplot(2,1,2);
    for i=1:length(SummaryCell_plot2);
        yvals=SummaryCell_plot2{i};
        ymean=mean(yvals);
        ysem=lt_sem(yvals);
        
        lt_plot_bar(i, ymean, {'Errors', ysem});
    end
    
    lt_subtitle('single targ has duplicate for both plots');
    
    
    
end
%% ===== PLOT ALL ONTO ONE PLOT [BOTH SYLS] [NO DUPLICATE ACROSS PLOTS]

lt_figure; hold on;
plot_text=0;

% ++++ INITIATE SUBPLOTS

% ==== plot 1: singe targ --> diff dir --> same dir ---> diff dir
lt_subplot(3,1,1); hold on;
title('single --> incongr --> congr --> incongr');
xlim([0 5]); ylim([-0.2 1.2])
% ==== plot 2
lt_subplot(3,1,2); hold on;
title('single --> congr --> incongr -->');
xlim([0 5]); ylim([-0.2 1.2])

% ==== plot 2
lt_subplot(3,1,3); hold on;
title('single only');
xlim([0 5]); ylim([-0.2 1.2])

SummaryCell_plot1=cell(1, 4);
SummaryCell_plot2=cell(1, 4);
SummaryCell_plot3=cell(1, 4);


for i=1:NumTotalExpts
    
    birdname=GLOB_STRUCT.info(i).birdname;
    exptname=GLOB_STRUCT.info(i).exptname;
    
    %     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    %     title([birdname '-' exptname]);
    %     xlim([0 50]);
    %     ylim([-0.3 1.2]);
    
    % ==== COLLECT data
    TypeAll=[]; % 0 = single; 1 = bidir; 2 = samedir
    ConsolAll=[];
    stardatesAll=[];
    
    % ======= singleTarg?
    type=0;
    if ~isempty(GLOB_STRUCT.data.singleTarg(i).consolidation)
        
        consol=GLOB_STRUCT.data.singleTarg(i).consolidation;
        startday=GLOB_STRUCT.data.singleTarg(i).consolidation_day1_relWN;
        
        % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
        
    end
    
    % ======= bidir?
    type = 1;
    if ~isempty(GLOB_STRUCT.data.multiDir(i).consolidation)
        
        consol=GLOB_STRUCT.data.multiDir(i).consolidation;
        startday=GLOB_STRUCT.data.multiDir(i).phase_day1_relWN;
        
        % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
    end
    
    % ======= bidir (second phase)?
    type = 1;
    if isfield(GLOB_STRUCT.data, 'multiDir2')
        
        if ~isempty(GLOB_STRUCT.data.multiDir2(i).consolidation)
            
            consol=GLOB_STRUCT.data.multiDir2(i).consolidation;
            startday=GLOB_STRUCT.data.multiDir2(i).phase_day1_relWN;
            
            % -- collect
            TypeAll=[TypeAll type];
            ConsolAll=[ConsolAll consol];
            stardatesAll=[stardatesAll startday];
        end
    end
    
    % ======= samedir?
    type = 2;
    if ~isempty(GLOB_STRUCT.data.sameDir(i).consolidation)
        
        consol=GLOB_STRUCT.data.sameDir(i).consolidation;
        startday=GLOB_STRUCT.data.sameDir(i).phase_day1_relWN;
        
        % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
        
        
    end
    
    % ++++ sort by start date of epochs
    [~, inds]=sort(stardatesAll);
    
    TypeAll=TypeAll(inds);
    ConsolAll=ConsolAll(inds);
    stardatesAll=stardatesAll(inds);
    
    %     if strcmp(birdname, 'bk34bk68') & strcmp(exptname, 'SeqDepPitchLMAN3');
    %         keyboard
    %     end
    
    % ++++ determine which figure to plot in based on order of types
    plotsMade=0;
    
    % --- figure 1?
    tmp=num2str(TypeAll);
    typeallStr=regexprep(tmp,'[^\w'']','');
    
    ind=regexp('0121', typeallStr);
    if length(typeallStr)==1,
        % so does not match multiple phases
        ind=ind(1);
    end
    
    if ~isempty(ind) & ind<3;
        plotX=ind:ind+length(typeallStr)-1;
        
        if ~(length(plotX)==1 & plotX==1);
            subplot(3,1,1);
            
            plot(plotX, ConsolAll, 'o-');
            plotsMade=plotsMade+1;
            if plot_text==1;
                lt_plot_text(max(plotX), ConsolAll(end), [birdname(1:4) '-' exptname]);
            end
            
            % ---- COLLECT FOR MEANS
            for k=1:length(plotX);
                SummaryCell_plot1{plotX(k)}=[SummaryCell_plot1{plotX(k)} ConsolAll(k)];
            end
        end
    end
    
    % --- figure 2?
    tmp=num2str(TypeAll);
    typeallStr=regexprep(tmp,'[^\w'']','');
    
    ind=regexp('0212', typeallStr);
    if length(typeallStr)==1,
        % so does not match multiple phases
        ind=ind(1);
    end
    
    if ~isempty(ind) & ind<3;
        plotX=ind:ind+length(typeallStr)-1;
        if ~(length(plotX)==1 & plotX==1);
            
            subplot(3,1,2);
            
            plot(plotX, ConsolAll, 'o-');
            plotsMade=plotsMade+1;
            if plot_text==1;
                lt_plot_text(max(plotX), ConsolAll(end), [birdname(1:4) '-' exptname]);
            end
            
            % ---- COLLECT FOR MEANS
            for k=1:length(plotX);
                SummaryCell_plot2{plotX(k)}=[SummaryCell_plot2{plotX(k)} ConsolAll(k)];
            end
        end
    end
    
    
    % --- figure 3?
    tmp=num2str(TypeAll);
    typeallStr=regexprep(tmp,'[^\w'']','');
    
    ind=regexp('0', typeallStr);
    %     if length(typeallStr)==1,
    %         % so does not match multiple phases
    %         ind=ind(1);
    %     end
    
    if ~isempty(ind)
        disp(TypeAll)
        plotX=ind;
        subplot(3,1,3);
        
        plot(plotX, ConsolAll, 'o-');
        plotsMade=plotsMade+1;
        if plot_text==1;
            lt_plot_text(max(plotX), ConsolAll(end), [birdname(1:4) '-' exptname]);
        end
        
        % ---- COLLECT FOR MEANS
        for k=1:length(plotX);
            SummaryCell_plot3{plotX(k)}=[SummaryCell_plot3{plotX(k)} ConsolAll(k)];
        end
    end
    
    
    %     % ++++++++++++++++++
    %     [tmp, plotX]=ismember(TypeAll, [0 1 2 1]);
    %     if ~isempty(plotX) & issorted(plotX) & ~any(diff(plotX)-1) & plotX(1)<3
    %         %     [~, ~, plotX]=intersect(TypeAll, [0 1 2 1]);
    %         %     if ~isempty(plotX) & issorted(plotX) & ~any(diff(plotX)-1) & plotX(1)<3
    %         % means this expt follows that order. <3 is to make sure matches
    %         % the start of the cycle (so not redundant across plots). issorted
    %         % is making sure is in chronological order. any...diff.. is to make
    %         % sure they are adjacent (e.g. 1 2 3, not 1 3)
    %
    %         subplot(2,1,1);
    %
    %         plot(plotX, ConsolAll, 'o-');
    %         plotsMade=plotsMade+1;
    %         if plot_text==1;
    %             lt_plot_text(max(plotX), ConsolAll(end), [birdname(1:4) '-' exptname]);
    %         end
    %
    %         % ---- COLLECT FOR MEANS
    %         for k=1:length(plotX);
    %             SummaryCell_plot1{plotX(k)}=[SummaryCell_plot1{plotX(k)} ConsolAll(k)];
    %         end
    %
    %     end
    %
    %     % --- figure 2?
    %     [tmp, plotX]=ismember(TypeAll, [0 2 1 2]);
    %     if ~isempty(plotX) & issorted(plotX) & ~any(diff(plotX)-1) & plotX(1)<3
    %         %     [~, ~, plotX]=intersect(TypeAll, [0 2 1 2]);
    %         %     if ~isempty(plotX) & issorted(plotX) & ~any(diff(plotX)-1) & plotX(1)<3
    %         % means this expt follows that order. <3 is to make sure matches
    %         % the start of the cycle (so not redundant across plots)
    %
    %         subplot(2,1,2);
    %
    %         plot(plotX, ConsolAll, 'o-');
    %         plotsMade=plotsMade+1;
    %
    %         if plot_text==1;
    %             lt_plot_text(max(plotX), ConsolAll(end), [birdname(1:4) '-' exptname]);
    %         end
    %
    %
    %         % ---- COLLECT FOR MEANS
    %         for k=1:length(plotX);
    %             SummaryCell_plot2{plotX(k)}=[SummaryCell_plot2{plotX(k)} ConsolAll(k)];
    %         end
    %
    %     end
    % +++++++++++++++++++++++++++++++++++++++++++
    
    if plotsMade>1 & TypeAll~=0
        disp(['PROBLEM ' birdname '-' exptname ' - plotted >1 time, and was not simply single-targ plots']);
    end
    
    if plotsMade<1
        disp(['PROBLEM!! - ' birdname '-' exptname ' failed to find correct subplot']);
    end
    
end


% ===== OVERLAY BARS FOR MEANS
subplot(3,1,1);
for i=1:length(SummaryCell_plot1);
    yvals=SummaryCell_plot1{i};
    ymean=mean(yvals);
    ysem=lt_sem(yvals);
    
    lt_plot_bar(i, ymean, {'Errors', ysem});
end

subplot(3,1,2);
for i=1:length(SummaryCell_plot2);
    yvals=SummaryCell_plot2{i};
    ymean=mean(yvals);
    ysem=lt_sem(yvals);
    
    lt_plot_bar(i, ymean, {'Errors', ysem});
end

subplot(3,1,3);
for i=1:length(SummaryCell_plot3);
    yvals=SummaryCell_plot3{i};
    ymean=mean(yvals);
    ysem=lt_sem(yvals);
    
    lt_plot_bar(i, ymean, {'Errors', ysem});
end



lt_subtitle('both targs combined (consol)');


%% ===== PLOT ALL ONTO ONE PLOT [FIRST TARG]
if (0)
    lt_figure; hold on;
    plot_text=0;
    
    % ++++ INITIATE SUBPLOTS
    
    % ==== plot 1: singe targ --> diff dir --> same dir ---> diff dir
    lt_subplot(2,1,1); hold on;
    title('single --> incongr --> congr --> incongr');
    xlim([0 5]); ylim([-0.2 1.2])
    % ==== plot 2
    lt_subplot(2,1,2); hold on;
    title('single --> congr --> incongr -->');
    xlim([0 5]); ylim([-0.2 1.2])
    ylabel('first targ only');
    
    SummaryCell_plot1=cell(1, 4);
    SummaryCell_plot2=cell(1, 4);
    
    
    for i=1:NumTotalExpts
        
        birdname=GLOB_STRUCT.info(i).birdname;
        exptname=GLOB_STRUCT.info(i).exptname;
        
        %     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        %     title([birdname '-' exptname]);
        %     xlim([0 50]);
        %     ylim([-0.3 1.2]);
        
        % ==== COLLECT data
        TypeAll=[]; % 0 = single; 1 = bidir; 2 = samedir
        ConsolAll=[];
        stardatesAll=[];
        
        % ======= singleTarg?
        type=0;
        if ~isempty(GLOB_STRUCT.data.singleTarg(i).consolidation)
            
            consol=GLOB_STRUCT.data.singleTarg(i).consolidation;
            startday=GLOB_STRUCT.data.singleTarg(i).consolidation_day1_relWN;
            
            % -- collect
            TypeAll=[TypeAll type];
            ConsolAll=[ConsolAll consol];
            stardatesAll=[stardatesAll startday];
            
        end
        
        % ======= bidir?
        type = 1;
        if ~isempty(GLOB_STRUCT.data.multiDir(i).consolidation)
            
            consol=GLOB_STRUCT.data.multiDir(i).consol_first;
            startday=GLOB_STRUCT.data.multiDir(i).phase_day1_relWN;
            
            % -- collect
            TypeAll=[TypeAll type];
            ConsolAll=[ConsolAll consol];
            stardatesAll=[stardatesAll startday];
        end
        
        % ======= bidir (second phase)?
        type = 1;
        if ~isempty(GLOB_STRUCT.data.multiDir2(i).consolidation)
            
            consol=GLOB_STRUCT.data.multiDir2(i).consol_first;
            startday=GLOB_STRUCT.data.multiDir2(i).phase_day1_relWN;
            
            % -- collect
            TypeAll=[TypeAll type];
            ConsolAll=[ConsolAll consol];
            stardatesAll=[stardatesAll startday];
        end
        
        % ======= samedir?
        type = 2;
        if ~isempty(GLOB_STRUCT.data.sameDir(i).consolidation)
            
            consol=GLOB_STRUCT.data.sameDir(i).consol_first;
            startday=GLOB_STRUCT.data.sameDir(i).phase_day1_relWN;
            
            % -- collect
            TypeAll=[TypeAll type];
            ConsolAll=[ConsolAll consol];
            stardatesAll=[stardatesAll startday];
            
            
        end
        
        % ++++ sort by start date of epochs
        [~, inds]=sort(stardatesAll);
        
        TypeAll=TypeAll(inds);
        ConsolAll=ConsolAll(inds);
        stardatesAll=stardatesAll(inds);
        
        %     if strcmp(birdname, 'bk34bk68') & strcmp(exptname, 'SeqDepPitchLMAN3');
        %         keyboard
        %     end
        
        % ++++ determine which figure to plot in based on order of types
        plotsMade=0;
        
        % --- figure 1?
        tmp=num2str(TypeAll);
        typeallStr=regexprep(tmp,'[^\w'']','');
        
        ind=regexp('0121', typeallStr);
        if length(typeallStr)==1,
            % so does not match multiple phases
            ind=ind(1);
        end
        
        if ~isempty(ind) & ind<3;
            plotX=ind:ind+length(typeallStr)-1;
            
            subplot(2,1,1);
            
            plot(plotX, ConsolAll, 'o-');
            plotsMade=plotsMade+1;
            if plot_text==1;
                lt_plot_text(max(plotX), ConsolAll(end), [birdname(1:4) '-' exptname]);
            end
            
            % ---- COLLECT FOR MEANS
            for k=1:length(plotX);
                SummaryCell_plot1{plotX(k)}=[SummaryCell_plot1{plotX(k)} ConsolAll(k)];
            end
            
        end
        
        % --- figure 2?
        tmp=num2str(TypeAll);
        typeallStr=regexprep(tmp,'[^\w'']','');
        
        ind=regexp('0212', typeallStr);
        if length(typeallStr)==1,
            % so does not match multiple phases
            ind=ind(1);
        end
        
        if ~isempty(ind) & ind<3;
            plotX=ind:ind+length(typeallStr)-1;
            
            subplot(2,1,2);
            
            plot(plotX, ConsolAll, 'o-');
            plotsMade=plotsMade+1;
            if plot_text==1;
                lt_plot_text(max(plotX), ConsolAll(end), [birdname(1:4) '-' exptname]);
            end
            
            % ---- COLLECT FOR MEANS
            for k=1:length(plotX);
                SummaryCell_plot2{plotX(k)}=[SummaryCell_plot2{plotX(k)} ConsolAll(k)];
            end
            
        end
        
        
        %     % ++++++++++++++++++
        %     [tmp, plotX]=ismember(TypeAll, [0 1 2 1]);
        %     if ~isempty(plotX) & issorted(plotX) & ~any(diff(plotX)-1) & plotX(1)<3
        %         %     [~, ~, plotX]=intersect(TypeAll, [0 1 2 1]);
        %         %     if ~isempty(plotX) & issorted(plotX) & ~any(diff(plotX)-1) & plotX(1)<3
        %         % means this expt follows that order. <3 is to make sure matches
        %         % the start of the cycle (so not redundant across plots). issorted
        %         % is making sure is in chronological order. any...diff.. is to make
        %         % sure they are adjacent (e.g. 1 2 3, not 1 3)
        %
        %         subplot(2,1,1);
        %
        %         plot(plotX, ConsolAll, 'o-');
        %         plotsMade=plotsMade+1;
        %         if plot_text==1;
        %             lt_plot_text(max(plotX), ConsolAll(end), [birdname(1:4) '-' exptname]);
        %         end
        %
        %         % ---- COLLECT FOR MEANS
        %         for k=1:length(plotX);
        %             SummaryCell_plot1{plotX(k)}=[SummaryCell_plot1{plotX(k)} ConsolAll(k)];
        %         end
        %
        %     end
        %
        %     % --- figure 2?
        %     [tmp, plotX]=ismember(TypeAll, [0 2 1 2]);
        %     if ~isempty(plotX) & issorted(plotX) & ~any(diff(plotX)-1) & plotX(1)<3
        %         %     [~, ~, plotX]=intersect(TypeAll, [0 2 1 2]);
        %         %     if ~isempty(plotX) & issorted(plotX) & ~any(diff(plotX)-1) & plotX(1)<3
        %         % means this expt follows that order. <3 is to make sure matches
        %         % the start of the cycle (so not redundant across plots)
        %
        %         subplot(2,1,2);
        %
        %         plot(plotX, ConsolAll, 'o-');
        %         plotsMade=plotsMade+1;
        %
        %         if plot_text==1;
        %             lt_plot_text(max(plotX), ConsolAll(end), [birdname(1:4) '-' exptname]);
        %         end
        %
        %
        %         % ---- COLLECT FOR MEANS
        %         for k=1:length(plotX);
        %             SummaryCell_plot2{plotX(k)}=[SummaryCell_plot2{plotX(k)} ConsolAll(k)];
        %         end
        %
        %     end
        % +++++++++++++++++++++++++++++++++++++++++++
        
        if plotsMade>1 & TypeAll~=0
            disp('PROBLEM - plotted >1 time, and was not simply single-targ plots');
        end
        
        if plotsMade<1
            disp(['PROBLEM!! - ' birdname '-' exptname ' failed to find correct subplot']);
        end
        
    end
    
    
    % ===== OVERLAY BARS FOR MEANS
    subplot(2,1,1);
    for i=1:length(SummaryCell_plot1);
        yvals=SummaryCell_plot1{i};
        ymean=mean(yvals);
        ysem=lt_sem(yvals);
        
        lt_plot_bar(i, ymean, {'Errors', ysem});
    end
    
    subplot(2,1,2);
    for i=1:length(SummaryCell_plot2);
        yvals=SummaryCell_plot2{i};
        ymean=mean(yvals);
        ysem=lt_sem(yvals);
        
        lt_plot_bar(i, ymean, {'Errors', ysem});
    end
    
    lt_subtitle('single targ has duplicate for both plots');
    
end
%% ===== PLOT ALL ONTO ONE PLOT [FIRST TARG SYLS] [NO DUPLICATE ACROSS PLOTS]

lt_figure; hold on;
plot_text=0;

% ++++ INITIATE SUBPLOTS

% ==== plot 1: singe targ --> diff dir --> same dir ---> diff dir
lt_subplot(3,1,1); hold on;
title('single --> incongr --> congr --> incongr');
xlim([0 5]); ylim([-0.2 1.2])
% ==== plot 2
lt_subplot(3,1,2); hold on;
title('single --> congr --> incongr -->');
xlim([0 5]); ylim([-0.2 1.2])

% ==== plot 2
lt_subplot(3,1,3); hold on;
title('single only');
xlim([0 5]); ylim([-0.2 1.2])

SummaryCell_plot1=cell(1, 4);
SummaryCell_plot2=cell(1, 4);
SummaryCell_plot3=cell(1, 4);


for i=1:NumTotalExpts
    
    birdname=GLOB_STRUCT.info(i).birdname;
    exptname=GLOB_STRUCT.info(i).exptname;
    
    %     [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    %     title([birdname '-' exptname]);
    %     xlim([0 50]);
    %     ylim([-0.3 1.2]);
    
    % ==== COLLECT data
    TypeAll=[]; % 0 = single; 1 = bidir; 2 = samedir
    ConsolAll=[];
    stardatesAll=[];
    
    % ======= singleTarg?
    type=0;
    if ~isempty(GLOB_STRUCT.data.singleTarg(i).consolidation)
        
        consol=GLOB_STRUCT.data.singleTarg(i).consolidation;
        startday=GLOB_STRUCT.data.singleTarg(i).consolidation_day1_relWN;
        
        % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
        
    end
    
    % ======= bidir?
    type = 1;
    if ~isempty(GLOB_STRUCT.data.multiDir(i).consolidation)
        
        consol=GLOB_STRUCT.data.multiDir(i).consol_first;
        startday=GLOB_STRUCT.data.multiDir(i).phase_day1_relWN;
        
        % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
    end
    
    % ======= bidir (second phase)?
    type = 1;
    if isfield(GLOB_STRUCT.data, 'multiDir2')
        
        if ~isempty(GLOB_STRUCT.data.multiDir2(i).consolidation)
            
            consol=GLOB_STRUCT.data.multiDir2(i).consol_first;
            startday=GLOB_STRUCT.data.multiDir2(i).phase_day1_relWN;
            
            % -- collect
            TypeAll=[TypeAll type];
            ConsolAll=[ConsolAll consol];
            stardatesAll=[stardatesAll startday];
        end
    end
    
    % ======= samedir?
    type = 2;
    if ~isempty(GLOB_STRUCT.data.sameDir(i).consolidation)
        
        consol=GLOB_STRUCT.data.sameDir(i).consol_first;
        startday=GLOB_STRUCT.data.sameDir(i).phase_day1_relWN;
        
        % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
        
        
    end
    
    % ++++ sort by start date of epochs
    [~, inds]=sort(stardatesAll);
    
    TypeAll=TypeAll(inds);
    ConsolAll=ConsolAll(inds);
    stardatesAll=stardatesAll(inds);
    
    %     if strcmp(birdname, 'bk34bk68') & strcmp(exptname, 'SeqDepPitchLMAN3');
    %         keyboard
    %     end
    
    % ++++ determine which figure to plot in based on order of types
    plotsMade=0;
    
    % --- figure 1?
    tmp=num2str(TypeAll);
    typeallStr=regexprep(tmp,'[^\w'']','');
    
    ind=regexp('0121', typeallStr);
    if length(typeallStr)==1,
        % so does not match multiple phases
        ind=ind(1);
    end
    
    if ~isempty(ind) & ind<3;
        plotX=ind:ind+length(typeallStr)-1;
        
        if ~(length(plotX)==1 & plotX==1);
            subplot(3,1,1);
            
            plot(plotX, ConsolAll, 'o-');
            plotsMade=plotsMade+1;
            if plot_text==1;
                lt_plot_text(max(plotX), ConsolAll(end), [birdname(1:4) '-' exptname]);
            end
            
            % ---- COLLECT FOR MEANS
            for k=1:length(plotX);
                SummaryCell_plot1{plotX(k)}=[SummaryCell_plot1{plotX(k)} ConsolAll(k)];
            end
        end
    end
    
    
    % --- figure 2?
    tmp=num2str(TypeAll);
    typeallStr=regexprep(tmp,'[^\w'']','');
    
    ind=regexp('0212', typeallStr);
    if length(typeallStr)==1,
        % so does not match multiple phases
        ind=ind(1);
    end
    
    if ~isempty(ind) & ind<3;
        plotX=ind:ind+length(typeallStr)-1;
        if ~(length(plotX)==1 & plotX==1);
            
            subplot(3,1,2);
            
            plot(plotX, ConsolAll, 'o-');
            plotsMade=plotsMade+1;
            if plot_text==1;
                lt_plot_text(max(plotX), ConsolAll(end), [birdname(1:4) '-' exptname]);
            end
            
            % ---- COLLECT FOR MEANS
            for k=1:length(plotX);
                SummaryCell_plot2{plotX(k)}=[SummaryCell_plot2{plotX(k)} ConsolAll(k)];
            end
        end
    end
    
    
    % --- figure 3?
    tmp=num2str(TypeAll);
    typeallStr=regexprep(tmp,'[^\w'']','');
    
    ind=regexp('0', typeallStr);
    %     if length(typeallStr)==1,
    %         % so does not match multiple phases
    %         ind=ind(1);
    %     end
    
    if ~isempty(ind)
        disp(TypeAll)
        plotX=ind;
        subplot(3,1,3);
        
        plot(plotX, ConsolAll, 'o-');
        plotsMade=plotsMade+1;
        if plot_text==1;
            lt_plot_text(max(plotX), ConsolAll(end), [birdname(1:4) '-' exptname]);
        end
        
        % ---- COLLECT FOR MEANS
        for k=1:length(plotX);
            SummaryCell_plot3{plotX(k)}=[SummaryCell_plot3{plotX(k)} ConsolAll(k)];
        end
    end
    
    
    %     % ++++++++++++++++++
    %     [tmp, plotX]=ismember(TypeAll, [0 1 2 1]);
    %     if ~isempty(plotX) & issorted(plotX) & ~any(diff(plotX)-1) & plotX(1)<3
    %         %     [~, ~, plotX]=intersect(TypeAll, [0 1 2 1]);
    %         %     if ~isempty(plotX) & issorted(plotX) & ~any(diff(plotX)-1) & plotX(1)<3
    %         % means this expt follows that order. <3 is to make sure matches
    %         % the start of the cycle (so not redundant across plots). issorted
    %         % is making sure is in chronological order. any...diff.. is to make
    %         % sure they are adjacent (e.g. 1 2 3, not 1 3)
    %
    %         subplot(2,1,1);
    %
    %         plot(plotX, ConsolAll, 'o-');
    %         plotsMade=plotsMade+1;
    %         if plot_text==1;
    %             lt_plot_text(max(plotX), ConsolAll(end), [birdname(1:4) '-' exptname]);
    %         end
    %
    %         % ---- COLLECT FOR MEANS
    %         for k=1:length(plotX);
    %             SummaryCell_plot1{plotX(k)}=[SummaryCell_plot1{plotX(k)} ConsolAll(k)];
    %         end
    %
    %     end
    %
    %     % --- figure 2?
    %     [tmp, plotX]=ismember(TypeAll, [0 2 1 2]);
    %     if ~isempty(plotX) & issorted(plotX) & ~any(diff(plotX)-1) & plotX(1)<3
    %         %     [~, ~, plotX]=intersect(TypeAll, [0 2 1 2]);
    %         %     if ~isempty(plotX) & issorted(plotX) & ~any(diff(plotX)-1) & plotX(1)<3
    %         % means this expt follows that order. <3 is to make sure matches
    %         % the start of the cycle (so not redundant across plots)
    %
    %         subplot(2,1,2);
    %
    %         plot(plotX, ConsolAll, 'o-');
    %         plotsMade=plotsMade+1;
    %
    %         if plot_text==1;
    %             lt_plot_text(max(plotX), ConsolAll(end), [birdname(1:4) '-' exptname]);
    %         end
    %
    %
    %         % ---- COLLECT FOR MEANS
    %         for k=1:length(plotX);
    %             SummaryCell_plot2{plotX(k)}=[SummaryCell_plot2{plotX(k)} ConsolAll(k)];
    %         end
    %
    %     end
    % +++++++++++++++++++++++++++++++++++++++++++
    
    if plotsMade>1 & TypeAll~=0
        disp(['PROBLEM ' birdname '-' exptname ' - plotted >1 time, and was not simply single-targ plots']);
    end
    
    if plotsMade<1
        disp(['PROBLEM!! - ' birdname '-' exptname ' failed to find correct subplot']);
    end
    
end


% ===== OVERLAY BARS FOR MEANS
subplot(3,1,1);
for i=1:length(SummaryCell_plot1);
    yvals=SummaryCell_plot1{i};
    ymean=mean(yvals);
    ysem=lt_sem(yvals);
    
    lt_plot_bar(i, ymean, {'Errors', ysem});
end

subplot(3,1,2);
for i=1:length(SummaryCell_plot2);
    yvals=SummaryCell_plot2{i};
    ymean=mean(yvals);
    ysem=lt_sem(yvals);
    
    lt_plot_bar(i, ymean, {'Errors', ysem});
end

subplot(3,1,3);
for i=1:length(SummaryCell_plot3);
    yvals=SummaryCell_plot3{i};
    ymean=mean(yvals);
    ysem=lt_sem(yvals);
    
    lt_plot_bar(i, ymean, {'Errors', ysem});
end


lt_subtitle('first targ only (consol)');


% ==============================  ANOVA - effect of group (DEP vs. IND)
% regardless of phase num?
% ----- NOTE THIS CODE NOT GENERAL, JUST USE FOR THIS PART

Yall=[];
CondAll=[];
PhaseAll=[];

% --- DEP, phase 2
vals=SummaryCell_plot1{2}; 
cond=1; % 1 = DEP; 2=IND
phase=2;
% -
Yall=[Yall vals];
CondAll=[CondAll cond*ones(1, length(vals))];
PhaseAll=[PhaseAll phase*ones(1,length(vals))];


% --- DEP, phase 3
vals=SummaryCell_plot2{3}; 
cond=1;
phase=3;
% -
Yall=[Yall vals];
CondAll=[CondAll cond*ones(1, length(vals))];
PhaseAll=[PhaseAll phase*ones(1,length(vals))];


% --- IND, phase 2
vals=SummaryCell_plot2{2}; 
cond=2;
phase=2;
% -
Yall=[Yall vals];
CondAll=[CondAll cond*ones(1, length(vals))];
PhaseAll=[PhaseAll phase*ones(1,length(vals))];


% ---- IND, phase 3
vals=SummaryCell_plot1{3}; 
cond=2;
phase=3;
% -
Yall=[Yall vals];
CondAll=[CondAll cond*ones(1, length(vals))];
PhaseAll=[PhaseAll phase*ones(1,length(vals))];

% --- RUN ANOVA
Group={};
Group{1}=CondAll;
Group{2}=PhaseAll;
anovan(Yall, Group);

inds=find(CondAll==1);
depvals=Yall(inds);

inds=find(CondAll==2);
indvals=Yall(inds);

lt_figure; hold on;
plot(1, depvals, 'ok');
plot(2, indvals, 'ok');
xlim([0 3]);


%% ===== [BOTH TARG COMBINED] PLOT ALL ONTO ONE PLOT [COMBINE ACROSS ALL EXPERIMENTAL TYPES]


% ==== COLLECT data
TypeAll=[]; % 0 = single; 1 = bidir; 2 = samedir
ConsolAll=[];
stardatesAll=[];
BirdnamesAll={};


for i=1:NumTotalExpts
    
    birdname=GLOB_STRUCT.info(i).birdname;
    exptname=GLOB_STRUCT.info(i).exptname;
    
    
    
    % ======= singleTarg?
    type=0;
    if ~isempty(GLOB_STRUCT.data.singleTarg(i).consolidation)
        
        consol=GLOB_STRUCT.data.singleTarg(i).consolidation;
        startday=GLOB_STRUCT.data.singleTarg(i).consolidation_day1_relWN;
        birdname=GLOB_STRUCT.data.singleTarg(i).birdname;
        BirdnamesAll=[BirdnamesAll birdname];
        
% -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
    end
    
    % ======= bidir?
    type = 1;
    if ~isempty(GLOB_STRUCT.data.multiDir(i).consolidation)
        
        consol=GLOB_STRUCT.data.multiDir(i).consolidation;
        startday=GLOB_STRUCT.data.multiDir(i).phase_day1_relWN;
        birdname=GLOB_STRUCT.data.multiDir(i).birdname;
        BirdnamesAll=[BirdnamesAll birdname];
        
        % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
    end
    
    % ======= bidir (second phase)?
    type = 1;
    if isfield(GLOB_STRUCT.data, 'multiDir2')
        
        if ~isempty(GLOB_STRUCT.data.multiDir2(i).consolidation)
            
            consol=GLOB_STRUCT.data.multiDir2(i).consolidation;
            startday=GLOB_STRUCT.data.multiDir2(i).phase_day1_relWN;
            birdname=GLOB_STRUCT.data.multiDir2(i).birdname;
            BirdnamesAll=[BirdnamesAll birdname];
           
            % -- collect
            TypeAll=[TypeAll type];
            ConsolAll=[ConsolAll consol];
            stardatesAll=[stardatesAll startday];
        end
    end
    
    % ======= samedir?
    type = 2;
    if ~isempty(GLOB_STRUCT.data.sameDir(i).consolidation)
        
        consol=GLOB_STRUCT.data.sameDir(i).consolidation;
        startday=GLOB_STRUCT.data.sameDir(i).phase_day1_relWN;
        birdname=GLOB_STRUCT.data.sameDir(i).birdname;
        BirdnamesAll=[BirdnamesAll birdname];
        
        % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
        
        
    end
end


% ================= PLOT
lt_figure; hold on;
title('consolid - both targs combined');

Yall={};

% -- single targ
x=1;
type=0;
color='k';

inds=TypeAll==type;
y=ConsolAll(inds);

plot(x, y, 'o', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors',lt_sem(y)});
Yall{1}=y;

% -- multidir
x=2;
type=1;
color='r';

inds=TypeAll==type;
y=ConsolAll(inds);

plot(x, y, 'o', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors',lt_sem(y)});
Yall{2}=y;

% -- samedir
x=3;
type=2;
color='b';

inds=TypeAll==type;
y=ConsolAll(inds);

plot(x, y, 'o', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors',lt_sem(y)});
Yall{3}=y;


% ==== paired lines
BirdnamesUnique=unique(BirdnamesAll);
for bb=1:length(BirdnamesUnique)
    birdname=BirdnamesUnique{bb};
    for i=0:2;
        ii=i+1;
        ind1=find(strcmp(BirdnamesAll, birdname) & TypeAll==i);
        ind2=find(strcmp(BirdnamesAll, birdname) & TypeAll==ii);
        
        if ~isempty(ind1) & ~isempty(ind2);
            for k=1:length(ind1);
                for kk=1:length(ind2);
                    X=[i+1 ii+1];
                    Y=[ConsolAll(ind1(k)) ConsolAll(ind2(kk))];
                    
                    plot(X, Y, '--k');
                end
            end
        end
    end
end

% === stats
for i=1:length(Yall)
    
    for ii=i+1:length(Yall)
        
        p=ranksum(Yall{i}, Yall{ii});
        %         [~, p]=ttest2(Yall{i}, Yall{ii});
        
        disp([num2str(i) ' - ' num2str(ii) '; p=' num2str(p)]);
        
        if p<0.2
            lt_plot_text(i, 1.1*max(Yall{i}), [num2str(i) ' - ' num2str(ii) '; p=' num2str(p)], 'r');
        end
    end
end


%% ===== [FIRST TARG] PLOT ALL ONTO ONE PLOT [COMBINE ACROSS ALL EXPERIMENTAL TYPES]


% ==== COLLECT data
TypeAll=[]; % 0 = single; 1 = bidir; 2 = samedir
ConsolAll=[];
stardatesAll=[];
ExptID=[];
BirdnamesAll={};


for i=1:NumTotalExpts
    
    birdname=GLOB_STRUCT.info(i).birdname;
    exptname=GLOB_STRUCT.info(i).exptname;
    
    
    
    % ======= singleTarg?
    type=0;
    if ~isempty(GLOB_STRUCT.data.singleTarg(i).consolidation)
        
        consol=GLOB_STRUCT.data.singleTarg(i).consolidation;
        startday=GLOB_STRUCT.data.singleTarg(i).consolidation_day1_relWN;
        
        birdname=GLOB_STRUCT.data.singleTarg(i).birdname;
        BirdnamesAll=[BirdnamesAll birdname];
        % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
        ExptID=[ExptID i];
        
    end
    
    % ======= bidir?
    type = 1;
    if ~isempty(GLOB_STRUCT.data.multiDir(i).consolidation)
        
        consol=GLOB_STRUCT.data.multiDir(i).consol_first;
        startday=GLOB_STRUCT.data.multiDir(i).phase_day1_relWN;
                birdname=GLOB_STRUCT.data.multiDir(i).birdname;
        BirdnamesAll=[BirdnamesAll birdname];

        % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
        ExptID=[ExptID i];
        
    end
    
    % ======= bidir (second phase)?
    type = 1;
    if isfield(GLOB_STRUCT.data, 'multiDir2')
        
        if ~isempty(GLOB_STRUCT.data.multiDir2(i).consolidation)
            
            consol=GLOB_STRUCT.data.multiDir2(i).consol_first;
            startday=GLOB_STRUCT.data.multiDir2(i).phase_day1_relWN;
            
            birdname=GLOB_STRUCT.data.multiDir2(i).birdname;
            BirdnamesAll=[BirdnamesAll birdname];
            % -- collect
            TypeAll=[TypeAll type];
            ConsolAll=[ConsolAll consol];
            stardatesAll=[stardatesAll startday];
            ExptID=[ExptID i];
            
        end
    end
    
    
    % ======= samedir?
    type = 2;
    if ~isempty(GLOB_STRUCT.data.sameDir(i).consolidation)
        
        consol=GLOB_STRUCT.data.sameDir(i).consol_first;
        startday=GLOB_STRUCT.data.sameDir(i).phase_day1_relWN;
        
        birdname=GLOB_STRUCT.data.sameDir(i).birdname;
        BirdnamesAll=[BirdnamesAll birdname];
        % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
        ExptID=[ExptID i];
        
        
        
    end
end


% ================= PLOT
lt_figure; hold on;
title('[IMPORTANT] first target only');
ylabel('consolidation');
xlabel('single targ -- incongr -- congr');
Yall={};

% -- single targ
x=1;
type=0;
color='k';

inds=TypeAll==type;
y=ConsolAll(inds);

plot(x, y, 'o', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors',lt_sem(y)});
Yall{1}=y;
% lt_plot_text(x, 1.4*mean(y), ['n=' num2str(length(y)) '; mean(sem)=' num2str(mean(y)) '(' num2str(lt_sem(y)) ')'])
lt_plot_text(x, 1.5*mean(y), ['n=' num2str(length(y)) '; mean(sem)=' num2str(mean(y)) '(' num2str(lt_sem(y)) ')'])

% -- multidir
x=2;
type=1;
color='r';

inds=TypeAll==type;
y=ConsolAll(inds);

plot(x, y, 'o', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors',lt_sem(y)});
Yall{2}=y;
lt_plot_text(x, 1.4*mean(y), ['n=' num2str(length(y)) '; mean(sem)=' num2str(mean(y)) '(' num2str(lt_sem(y)) ')'])

% -- samedir
x=3;
type=2;
color='b';

inds=TypeAll==type;
y=ConsolAll(inds);

plot(x, y, 'o', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors',lt_sem(y)});
Yall{3}=y;
lt_plot_text(x, 1.4*mean(y), ['n=' num2str(length(y)) '; mean(sem)=' num2str(mean(y)) '(' num2str(lt_sem(y)) ')'])

% ==== paired lines
BirdnamesUnique=unique(BirdnamesAll);
for bb=1:length(BirdnamesUnique)
    birdname=BirdnamesUnique{bb};
    for i=0:2;
        ii=i+1;
        ind1=find(strcmp(BirdnamesAll, birdname) & TypeAll==i);
        ind2=find(strcmp(BirdnamesAll, birdname) & TypeAll==ii);
        
        if ~isempty(ind1) & ~isempty(ind2);
            for k=1:length(ind1);
                for kk=1:length(ind2);
                    X=[i+1 ii+1];
                    Y=[ConsolAll(ind1(k)) ConsolAll(ind2(kk))];
                    
                    plot(X, Y, '--k');
                end
            end
        end
    end
end




% === stats
for i=1:length(Yall)
    
    for ii=i+1:length(Yall)
        
        p=ranksum(Yall{i}, Yall{ii});
        %         [~, p]=ttest2(Yall{i}, Yall{ii});
        
        disp([num2str(i) ' - ' num2str(ii) '; p=' num2str(p)]);
        
        if p<0.2
            lt_plot_text(i, 1.1*max(Yall{i}), [num2str(i) ' - ' num2str(ii) '; (ranksum)p=' num2str(p)], 'r');
        end
    end
end
%% ===== [AFP BIAS ONLY] PLOT ALL ONTO ONE PLOT [COMBINE ACROSS ALL EXPERIMENTAL TYPES]


% ==== COLLECT data
TypeAll=[]; % 0 = single; 1 = bidir; 2 = samedir
ConsolAll=[];
stardatesAll=[];
BirdnamesAll={};


for i=1:NumTotalExpts
    
    birdname=GLOB_STRUCT.info(i).birdname;
    exptname=GLOB_STRUCT.info(i).exptname;
    
    
    
    % ======= singleTarg?
    type=0;
    if ~isempty(GLOB_STRUCT.data.singleTarg(i).consolidation)
        
        consol=GLOB_STRUCT.data.singleTarg(i).ffPBS-GLOB_STRUCT.data.singleTarg(i).ffMUSC;
        startday=GLOB_STRUCT.data.singleTarg(i).consolidation_day1_relWN;
                birdname=GLOB_STRUCT.data.singleTarg(i).birdname;
        BirdnamesAll=[BirdnamesAll birdname];

        % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
        
    end
    
    % ======= bidir?
    type = 1;
    if ~isempty(GLOB_STRUCT.data.multiDir(i).consolidation)
        
        consol=GLOB_STRUCT.data.multiDir(i).ff_combined_PBS-GLOB_STRUCT.data.multiDir(i).ff_combined_MUSC;
        startday=GLOB_STRUCT.data.multiDir(i).phase_day1_relWN;
                        birdname=GLOB_STRUCT.data.multiDir(i).birdname;
        BirdnamesAll=[BirdnamesAll birdname];

        % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
    end
    
    % ======= bidir (second phase)?
    type = 1;
    if isfield(GLOB_STRUCT.data, 'multiDir2')
        
        if ~isempty(GLOB_STRUCT.data.multiDir2(i).consolidation)
            
            consol=GLOB_STRUCT.data.multiDir2(i).ff_combined_PBS-GLOB_STRUCT.data.multiDir2(i).ff_combined_MUSC;
            startday=GLOB_STRUCT.data.multiDir2(i).phase_day1_relWN;
            
             birdname=GLOB_STRUCT.data.multiDir2(i).birdname;
            BirdnamesAll=[BirdnamesAll birdname];
           % -- collect
            TypeAll=[TypeAll type];
            ConsolAll=[ConsolAll consol];
            stardatesAll=[stardatesAll startday];
        end
    end
    
    % ======= samedir?
    type = 2;
    if ~isempty(GLOB_STRUCT.data.sameDir(i).consolidation)
        
        consol=GLOB_STRUCT.data.sameDir(i).ff_combined_PBS-GLOB_STRUCT.data.sameDir(i).ff_combined_MUSC;
        startday=GLOB_STRUCT.data.sameDir(i).phase_day1_relWN;
        
        birdname=GLOB_STRUCT.data.sameDir(i).birdname;
        BirdnamesAll=[BirdnamesAll birdname];
        % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
        
        
    end
end


% ================= PLOT
lt_figure; hold on;
title('AFP bias [both targ combine]');

Yall={};

% -- single targ
x=1;
type=0;
color='k';

inds=TypeAll==type;
y=ConsolAll(inds);

plot(x, y, 'o', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors',lt_sem(y)});
Yall{1}=y;

% -- multidir
x=2;
type=1;
color='r';

inds=TypeAll==type;
y=ConsolAll(inds);

plot(x, y, 'o', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors',lt_sem(y)});
Yall{2}=y;

% -- samedir
x=3;
type=2;
color='b';

inds=TypeAll==type;
y=ConsolAll(inds);

plot(x, y, 'o', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors',lt_sem(y)});
Yall{3}=y;

% ==== paired lines
BirdnamesUnique=unique(BirdnamesAll);
for bb=1:length(BirdnamesUnique)
    birdname=BirdnamesUnique{bb};
    for i=0:2;
        ii=i+1;
        ind1=find(strcmp(BirdnamesAll, birdname) & TypeAll==i);
        ind2=find(strcmp(BirdnamesAll, birdname) & TypeAll==ii);
        
        if ~isempty(ind1) & ~isempty(ind2);
            for k=1:length(ind1);
                for kk=1:length(ind2);
                    X=[i+1 ii+1];
                    Y=[ConsolAll(ind1(k)) ConsolAll(ind2(kk))];
                    
                    plot(X, Y, '--k');
                end
            end
        end
    end
end



% === stats
for i=1:length(Yall)
    
    for ii=i+1:length(Yall)
        
        p=ranksum(Yall{i}, Yall{ii});
        %         [~, p]=ttest2(Yall{i}, Yall{ii});
        
        disp([num2str(i) ' - ' num2str(ii) '; p=' num2str(p)]);
        
        if p<0.2
            lt_plot_text(i, 1.1*max(Yall{i}), [num2str(i) ' - ' num2str(ii) '; p=' num2str(p)], 'r');
        end
    end
end

%% ===== [FIRST TARG] [AFP BIAS ONLY] PLOT ALL ONTO ONE PLOT [COMBINE ACROSS ALL EXPERIMENTAL TYPES]


% ==== COLLECT data
TypeAll=[]; % 0 = single; 1 = bidir; 2 = samedir
ConsolAll=[];
stardatesAll=[];
BirdnamesAll={};

for i=1:NumTotalExpts
    
    birdname=GLOB_STRUCT.info(i).birdname;
    exptname=GLOB_STRUCT.info(i).exptname;
    
    
    
    % ======= singleTarg?
    type=0;
    if ~isempty(GLOB_STRUCT.data.singleTarg(i).consolidation)
        
        consol=GLOB_STRUCT.data.singleTarg(i).ffPBS-GLOB_STRUCT.data.singleTarg(i).ffMUSC;
        startday=GLOB_STRUCT.data.singleTarg(i).consolidation_day1_relWN;
        
        birdname=GLOB_STRUCT.data.singleTarg(i).birdname;
        BirdnamesAll=[BirdnamesAll birdname];
        % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
        
    end
    
    % ======= bidir?
    type = 1;
    if ~isempty(GLOB_STRUCT.data.multiDir(i).consolidation)
        
        consol=GLOB_STRUCT.data.multiDir(i).ff_first_PBS-GLOB_STRUCT.data.multiDir(i).ff_first_MUSC;
        startday=GLOB_STRUCT.data.multiDir(i).phase_day1_relWN;
        
                 birdname=GLOB_STRUCT.data.multiDir(i).birdname;
        BirdnamesAll=[BirdnamesAll birdname];
       % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
    end
    
    % ======= bidir (second phase)?
    type = 1;
    if isfield(GLOB_STRUCT.data, 'multiDir2')
        
        if ~isempty(GLOB_STRUCT.data.multiDir2(i).consolidation)
            
            consol=GLOB_STRUCT.data.multiDir2(i).ff_first_PBS-GLOB_STRUCT.data.multiDir2(i).ff_first_MUSC;
            startday=GLOB_STRUCT.data.multiDir2(i).phase_day1_relWN;
            
              birdname=GLOB_STRUCT.data.multiDir2(i).birdname;
            BirdnamesAll=[BirdnamesAll birdname];
          % -- collect
            TypeAll=[TypeAll type];
            ConsolAll=[ConsolAll consol];
            stardatesAll=[stardatesAll startday];
        end
    end
    
    % ======= samedir?
    type = 2;
    if ~isempty(GLOB_STRUCT.data.sameDir(i).consolidation)
        
        consol=GLOB_STRUCT.data.sameDir(i).ff_first_PBS-GLOB_STRUCT.data.sameDir(i).ff_first_MUSC;
        startday=GLOB_STRUCT.data.sameDir(i).phase_day1_relWN;
        
         birdname=GLOB_STRUCT.data.sameDir(i).birdname;
        BirdnamesAll=[BirdnamesAll birdname];
       % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
        
        
    end
end


% ================= PLOT
lt_figure; hold on;
title('AFP bias [single targ]');

Yall={};

% -- single targ
x=1;
type=0;
color='k';

inds=TypeAll==type;
y=ConsolAll(inds);

plot(x, y, 'o', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors',lt_sem(y)});
Yall{1}=y;

% -- multidir
x=2;
type=1;
color='r';

inds=TypeAll==type;
y=ConsolAll(inds);

plot(x, y, 'o', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors',lt_sem(y)});
Yall{2}=y;

% -- samedir
x=3;
type=2;
color='b';

inds=TypeAll==type;
y=ConsolAll(inds);

plot(x, y, 'o', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors',lt_sem(y)});
Yall{3}=y;

% ==== paired lines
BirdnamesUnique=unique(BirdnamesAll);
for bb=1:length(BirdnamesUnique)
    birdname=BirdnamesUnique{bb};
    for i=0:2;
        ii=i+1;
        ind1=find(strcmp(BirdnamesAll, birdname) & TypeAll==i);
        ind2=find(strcmp(BirdnamesAll, birdname) & TypeAll==ii);
        
        if ~isempty(ind1) & ~isempty(ind2);
            for k=1:length(ind1);
                for kk=1:length(ind2);
                    X=[i+1 ii+1];
                    Y=[ConsolAll(ind1(k)) ConsolAll(ind2(kk))];
                    
                    plot(X, Y, '--k');
                end
            end
        end
    end
end

% === stats
for i=1:length(Yall)
    
    for ii=i+1:length(Yall)
        
        p=ranksum(Yall{i}, Yall{ii});
        %         [~, p]=ttest2(Yall{i}, Yall{ii});
        
        disp([num2str(i) ' - ' num2str(ii) '; p=' num2str(p)]);
        
        if p<0.2
            lt_plot_text(i, 1.1*max(Yall{i}), [num2str(i) ' - ' num2str(ii) '; p=' num2str(p)], 'r');
        end
    end
end


%% ===== [SECOND TARG] [AFP BIAS ONLY] PLOT ALL ONTO ONE PLOT [COMBINE ACROSS ALL EXPERIMENTAL TYPES]


% ==== COLLECT data
TypeAll=[]; % 0 = single; 1 = bidir; 2 = samedir
ConsolAll=[];
stardatesAll=[];
BirdnamesAll={};


for i=1:NumTotalExpts
    
    birdname=GLOB_STRUCT.info(i).birdname;
    exptname=GLOB_STRUCT.info(i).exptname;
    
    
    
%     % ======= singleTarg?
%     type=0;
%     if ~isempty(GLOB_STRUCT.data.singleTarg(i).consolidation)
%         
%         consol=GLOB_STRUCT.data.singleTarg(i).ffPBS-GLOB_STRUCT.data.singleTarg(i).ffMUSC;
%         startday=GLOB_STRUCT.data.singleTarg(i).consolidation_day1_relWN;
%         
%         % -- collect
%         TypeAll=[TypeAll type];
%         ConsolAll=[ConsolAll consol];
%         stardatesAll=[stardatesAll startday];
%         
%     end
    
    % ======= bidir?
    type = 1;
    if ~isempty(GLOB_STRUCT.data.multiDir(i).consolidation)
        
        consol=GLOB_STRUCT.data.multiDir(i).ff_second_PBS-GLOB_STRUCT.data.multiDir(i).ff_second_MUSC;
        startday=GLOB_STRUCT.data.multiDir(i).phase_day1_relWN;
        
                 birdname=GLOB_STRUCT.data.multiDir(i).birdname;
        BirdnamesAll=[BirdnamesAll birdname];
       % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
    end
    
    % ======= bidir (second phase)?
    type = 1;
    if isfield(GLOB_STRUCT.data, 'multiDir2')
        
        if ~isempty(GLOB_STRUCT.data.multiDir2(i).consolidation)
            
            consol=GLOB_STRUCT.data.multiDir2(i).ff_second_PBS-GLOB_STRUCT.data.multiDir2(i).ff_second_MUSC;
            startday=GLOB_STRUCT.data.multiDir2(i).phase_day1_relWN;
             birdname=GLOB_STRUCT.data.multiDir2(i).birdname;
            BirdnamesAll=[BirdnamesAll birdname];
           
            % -- collect
            TypeAll=[TypeAll type];
            ConsolAll=[ConsolAll consol];
            stardatesAll=[stardatesAll startday];
        end
    end
    
    % ======= samedir?
    type = 2;
    if ~isempty(GLOB_STRUCT.data.sameDir(i).consolidation)
        
        consol=GLOB_STRUCT.data.sameDir(i).ff_second_PBS-GLOB_STRUCT.data.sameDir(i).ff_second_MUSC;
        startday=GLOB_STRUCT.data.sameDir(i).phase_day1_relWN;
        
        birdname=GLOB_STRUCT.data.sameDir(i).birdname;
        BirdnamesAll=[BirdnamesAll birdname];
        % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
        
        
    end
end


% ================= PLOT
lt_figure; hold on;
title('AFP bias [second targ]');

Yall={};

% % -- single targ
% x=1;
% type=0;
% color='k';
% 
% inds=TypeAll==type;
% y=ConsolAll(inds);
% 
% plot(x, y, 'o', 'Color',color);
% lt_plot_bar(x, mean(y), {'Errors',lt_sem(y)});
% Yall{1}=y;

% -- multidir
x=2;
type=1;
color='r';

inds=TypeAll==type;
y=ConsolAll(inds);

plot(x, y, 'o', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors',lt_sem(y)});
Yall{2}=y;

% -- samedir
x=3;
type=2;
color='b';

inds=TypeAll==type;
y=ConsolAll(inds);

plot(x, y, 'o', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors',lt_sem(y)});
Yall{3}=y;

% ==== paired lines
BirdnamesUnique=unique(BirdnamesAll);
for bb=1:length(BirdnamesUnique)
    birdname=BirdnamesUnique{bb};
    for i=0:2;
        ii=i+1;
        ind1=find(strcmp(BirdnamesAll, birdname) & TypeAll==i);
        ind2=find(strcmp(BirdnamesAll, birdname) & TypeAll==ii);
        
        if ~isempty(ind1) & ~isempty(ind2);
            for k=1:length(ind1);
                for kk=1:length(ind2);
                    X=[i+1 ii+1];
                    Y=[ConsolAll(ind1(k)) ConsolAll(ind2(kk))];
                    
                    plot(X, Y, '--k');
                end
            end
        end
    end
end


% === stats
        p=ranksum(Yall{2}, Yall{3});
        %         [~, p]=ttest2(Yall{i}, Yall{ii});
                
        if p<0.2
            lt_plot_text(2, 1.1*max(Yall{2}), ['p=' num2str(p)], 'r');
        end

%% ===== [HIT RATE, PBS] PLOT ALL ONTO ONE PLOT [COMBINE ACROSS ALL EXPERIMENTAL TYPES]


% ==== COLLECT data
TypeAll=[]; % 0 = single; 1 = bidir; 2 = samedir
ConsolAll=[];
stardatesAll=[];
BirdnamesAll={};

for i=1:NumTotalExpts
    
    birdname=GLOB_STRUCT.info(i).birdname;
    exptname=GLOB_STRUCT.info(i).exptname;
    
    
    
    % ======= singleTarg?
    type=0;
    if ~isempty(GLOB_STRUCT.data.singleTarg(i).hrate_PBS)
        
        consol=GLOB_STRUCT.data.singleTarg(i).hrate_PBS;
        startday=GLOB_STRUCT.data.singleTarg(i).consolidation_day1_relWN;
        
                birdname=GLOB_STRUCT.data.singleTarg(i).birdname;
        BirdnamesAll=[BirdnamesAll birdname];
        % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
        
    end
    
    % ======= bidir?
    type = 1;
    if ~isempty(GLOB_STRUCT.data.multiDir(i).hrate_combined_PBS)
        
        consol=GLOB_STRUCT.data.multiDir(i).hrate_combined_PBS;
        startday=GLOB_STRUCT.data.multiDir(i).phase_day1_relWN;
                        birdname=GLOB_STRUCT.data.multiDir(i).birdname;
        BirdnamesAll=[BirdnamesAll birdname];
        
        % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
    end
    
    %     % ======= bidir (second phase)?
    %     type = 1;
    %     if ~isempty(GLOB_STRUCT.data.multiDir2(i).)
    %
    %         consol=GLOB_STRUCT.data.multiDir2(i).consolidation;
    %         startday=GLOB_STRUCT.data.multiDir2(i).phase_day1_relWN;
    %
    %         % -- collect
    %         TypeAll=[TypeAll type];
    %         ConsolAll=[ConsolAll consol];
    %         stardatesAll=[stardatesAll startday];
    %     end
    
    % ======= samedir?
    type = 2;
    if ~isempty(GLOB_STRUCT.data.sameDir(i).hrate_combined_PBS)
        
        consol=GLOB_STRUCT.data.sameDir(i).hrate_combined_PBS;
        startday=GLOB_STRUCT.data.sameDir(i).phase_day1_relWN;
        
         birdname=GLOB_STRUCT.data.sameDir(i).birdname;
        BirdnamesAll=[BirdnamesAll birdname];
       % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
        
        
    end
end


% ================= PLOT
lt_figure; hold on;
title('hit rate, PBS [combined]');
Yall={};

% -- single targ
x=1;
type=0;
color='k';

inds=TypeAll==type;
y=ConsolAll(inds);

plot(x, y, 'o', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors',lt_sem(y)});
Yall{1}=y;

% -- multidir
x=2;
type=1;
color='r';

inds=TypeAll==type;
y=ConsolAll(inds);

plot(x, y, 'o', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors',lt_sem(y)});
Yall{2}=y;

% -- samedir
x=3;
type=2;
color='b';

inds=TypeAll==type;
y=ConsolAll(inds);

plot(x, y, 'o', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors',lt_sem(y)});
Yall{3}=y;

% ==== paired lines
BirdnamesUnique=unique(BirdnamesAll);
for bb=1:length(BirdnamesUnique)
    birdname=BirdnamesUnique{bb};
    for i=0:2;
        ii=i+1;
        ind1=find(strcmp(BirdnamesAll, birdname) & TypeAll==i);
        ind2=find(strcmp(BirdnamesAll, birdname) & TypeAll==ii);
        
        if ~isempty(ind1) & ~isempty(ind2);
            for k=1:length(ind1);
                for kk=1:length(ind2);
                    X=[i+1 ii+1];
                    Y=[ConsolAll(ind1(k)) ConsolAll(ind2(kk))];
                    
                    plot(X, Y, '--k');
                end
            end
        end
    end
end

% === stats
for i=1:length(Yall)
    
    for ii=i+1:length(Yall)
        
        p=ranksum(Yall{i}, Yall{ii});
        %         [~, p]=ttest2(Yall{i}, Yall{ii});
        
        disp([num2str(i) ' - ' num2str(ii) '; p=' num2str(p)]);

        if p<0.2
            lt_plot_text(i, 1.1*max(Yall{i}), [num2str(i) ' - ' num2str(ii) '; p=' num2str(p)], 'r');
        end
        
    end
end


%% ===== [HIT RATE, MUSC] PLOT ALL ONTO ONE PLOT [COMBINE ACROSS ALL EXPERIMENTAL TYPES]


% ==== COLLECT data
TypeAll=[]; % 0 = single; 1 = bidir; 2 = samedir
ConsolAll=[];
stardatesAll=[];
BirdnamesAll={};

for i=1:NumTotalExpts
    
    birdname=GLOB_STRUCT.info(i).birdname;
    exptname=GLOB_STRUCT.info(i).exptname;
    
    
    
    % ======= singleTarg?
    type=0;
    if ~isempty(GLOB_STRUCT.data.singleTarg(i).hrate_MUSC)
        
        consol=GLOB_STRUCT.data.singleTarg(i).hrate_MUSC;
        startday=GLOB_STRUCT.data.singleTarg(i).consolidation_day1_relWN;
        birdname=GLOB_STRUCT.data.singleTarg(i).birdname;
        BirdnamesAll=[BirdnamesAll birdname];
        
        % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
        
    end
    
    % ======= bidir?
    type = 1;
    if ~isempty(GLOB_STRUCT.data.multiDir(i).hrate_combined_MUSC)
        
        consol=GLOB_STRUCT.data.multiDir(i).hrate_combined_MUSC;
        startday=GLOB_STRUCT.data.multiDir(i).phase_day1_relWN;
        
                 birdname=GLOB_STRUCT.data.multiDir(i).birdname;
        BirdnamesAll=[BirdnamesAll birdname];
        % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
    end
    
    %     % ======= bidir (second phase)?
    %     type = 1;
    %     if ~isempty(GLOB_STRUCT.data.multiDir2(i).)
    %
    %         consol=GLOB_STRUCT.data.multiDir2(i).consolidation;
    %         startday=GLOB_STRUCT.data.multiDir2(i).phase_day1_relWN;
    %
    %         % -- collect
    %         TypeAll=[TypeAll type];
    %         ConsolAll=[ConsolAll consol];
    %         stardatesAll=[stardatesAll startday];
    %     end
    
    % ======= samedir?
    type = 2;
    if ~isempty(GLOB_STRUCT.data.sameDir(i).hrate_combined_MUSC)
        
        consol=GLOB_STRUCT.data.sameDir(i).hrate_combined_MUSC;
        startday=GLOB_STRUCT.data.sameDir(i).phase_day1_relWN;
        
         birdname=GLOB_STRUCT.data.sameDir(i).birdname;
        BirdnamesAll=[BirdnamesAll birdname];
        % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
        
        
    end
end


% ================= PLOT
lt_figure; hold on;
title('hit rate, PBS');
Yall={};

% -- single targ
x=1;
type=0;
color='k';

inds=TypeAll==type;
y=ConsolAll(inds);

plot(x, y, 'o', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors',lt_sem(y)});
Yall{1}=y;

% -- multidir
x=2;
type=1;
color='r';

inds=TypeAll==type;
y=ConsolAll(inds);

plot(x, y, 'o', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors',lt_sem(y)});
Yall{2}=y;

% -- samedir
x=3;
type=2;
color='b';

inds=TypeAll==type;
y=ConsolAll(inds);

plot(x, y, 'o', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors',lt_sem(y)});
Yall{3}=y;

% ==== paired lines
BirdnamesUnique=unique(BirdnamesAll);
for bb=1:length(BirdnamesUnique)
    birdname=BirdnamesUnique{bb};
    for i=0:2;
        ii=i+1;
        ind1=find(strcmp(BirdnamesAll, birdname) & TypeAll==i);
        ind2=find(strcmp(BirdnamesAll, birdname) & TypeAll==ii);
        
        if ~isempty(ind1) & ~isempty(ind2);
            for k=1:length(ind1);
                for kk=1:length(ind2);
                    X=[i+1 ii+1];
                    Y=[ConsolAll(ind1(k)) ConsolAll(ind2(kk))];
                    
                    plot(X, Y, '--k');
                end
            end
        end
    end
end

% === stats
% === stats
for i=1:length(Yall)
    
    for ii=i+1:length(Yall)
        
        p=ranksum(Yall{i}, Yall{ii});
        %         [~, p]=ttest2(Yall{i}, Yall{ii});
        
        disp([num2str(i) ' - ' num2str(ii) '; p=' num2str(p)]);

        if p<0.2
            lt_plot_text(i, 1.1*max(Yall{i}), [num2str(i) ' - ' num2str(ii) '; p=' num2str(p)], 'r');
        end
        
    end
end

%% ===== [HIT RATE, DIFFERENCE (PBS - MUSC)] PLOT ALL ONTO ONE PLOT [COMBINE ACROSS ALL EXPERIMENTAL TYPES]


% ==== COLLECT data
TypeAll=[]; % 0 = single; 1 = bidir; 2 = samedir
ConsolAll=[];
stardatesAll=[];
BirdnamesAll={};

for i=1:NumTotalExpts
    
    birdname=GLOB_STRUCT.info(i).birdname;
    exptname=GLOB_STRUCT.info(i).exptname;
    
    
    
    % ======= singleTarg?
    type=0;
    if ~isempty(GLOB_STRUCT.data.singleTarg(i).hrate_PBS)
        
        consol=GLOB_STRUCT.data.singleTarg(i).hrate_MUSC-GLOB_STRUCT.data.singleTarg(i).hrate_PBS;
        startday=GLOB_STRUCT.data.singleTarg(i).consolidation_day1_relWN;
        
                  birdname=GLOB_STRUCT.data.singleTarg(i).birdname;
        BirdnamesAll=[BirdnamesAll birdname];
       % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
        
    end
    
    % ======= bidir?
    type = 1;
    if ~isempty(GLOB_STRUCT.data.multiDir(i).hrate_combined_PBS)
        
        consol=GLOB_STRUCT.data.multiDir(i).hrate_combined_MUSC - GLOB_STRUCT.data.multiDir(i).hrate_combined_PBS;
        startday=GLOB_STRUCT.data.multiDir(i).phase_day1_relWN;
              birdname=GLOB_STRUCT.data.multiDir(i).birdname;
            BirdnamesAll=[BirdnamesAll birdname];
       
        % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
    end
    
    %     % ======= bidir (second phase)?
    %     type = 1;
    %     if ~isempty(GLOB_STRUCT.data.multiDir2(i).)
    %
    %         consol=GLOB_STRUCT.data.multiDir2(i).consolidation;
    %         startday=GLOB_STRUCT.data.multiDir2(i).phase_day1_relWN;
    %
    %         % -- collect
    %         TypeAll=[TypeAll type];
    %         ConsolAll=[ConsolAll consol];
    %         stardatesAll=[stardatesAll startday];
    %     end
    
    % ======= samedir?
    type = 2;
    if ~isempty(GLOB_STRUCT.data.sameDir(i).hrate_combined_PBS)
        
        consol=GLOB_STRUCT.data.sameDir(i).hrate_combined_MUSC - GLOB_STRUCT.data.sameDir(i).hrate_combined_PBS;
        startday=GLOB_STRUCT.data.sameDir(i).phase_day1_relWN;
        
              birdname=GLOB_STRUCT.data.sameDir(i).birdname;
            BirdnamesAll=[BirdnamesAll birdname];
       % -- collect
        TypeAll=[TypeAll type];
        ConsolAll=[ConsolAll consol];
        stardatesAll=[stardatesAll startday];
        
        
    end
end


% ================= PLOT
lt_figure; hold on;
title('hit rate (PBS -  MUSC)');
Yall={};

% -- single targ
x=1;
type=0;
color='k';

inds=TypeAll==type;
y=ConsolAll(inds);

plot(x, y, 'o', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors',lt_sem(y)});
Yall{1}=y;

% -- multidir
x=2;
type=1;
color='r';

inds=TypeAll==type;
y=ConsolAll(inds);

plot(x, y, 'o', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors',lt_sem(y)});
Yall{2}=y;

% -- samedir
x=3;
type=2;
color='b';

inds=TypeAll==type;
y=ConsolAll(inds);

plot(x, y, 'o', 'Color',color);
lt_plot_bar(x, mean(y), {'Errors',lt_sem(y)});
Yall{3}=y;
% ==== paired lines
BirdnamesUnique=unique(BirdnamesAll);
for bb=1:length(BirdnamesUnique)
    birdname=BirdnamesUnique{bb};
    for i=0:2;
        ii=i+1;
        ind1=find(strcmp(BirdnamesAll, birdname) & TypeAll==i);
        ind2=find(strcmp(BirdnamesAll, birdname) & TypeAll==ii);
        
        if ~isempty(ind1) & ~isempty(ind2);
            for k=1:length(ind1);
                for kk=1:length(ind2);
                    X=[i+1 ii+1];
                    Y=[ConsolAll(ind1(k)) ConsolAll(ind2(kk))];
                    
                    plot(X, Y, '--k');
                end
            end
        end
    end
end


% === stats
for i=1:length(Yall)
    
    for ii=i+1:length(Yall)
        
        p=ranksum(Yall{i}, Yall{ii});
        %         [~, p]=ttest2(Yall{i}, Yall{ii});
        
        disp([num2str(i) ' - ' num2str(ii) '; p=' num2str(p)]);
        
    end
end

%% ================================ MUSC SHIFTS DRAGGING EACH OTHER DOWN?
lt_figure; hold on;

% ======= MUSC

% == Bidir, MUSC
lt_subplot(4,2,1); hold on;
title('bidir, MUSC');
xlabel('first ctxt');
ylabel('second ctxt');
bidir_First_MUSC=[GLOB_STRUCT.data.multiDir.ff_first_MUSC];
bidir_Second_MUSC=[GLOB_STRUCT.data.multiDir.ff_second_MUSC];

lt_regress(bidir_Second_MUSC, bidir_First_MUSC, 1, 0, 1, 1, 'r');
lt_plot_zeroline;
lt_plot_zeroline_vert;

% == Bidir, AFP
lt_subplot(4,2,2); hold on;
title('bidir, AFP');

bidir_First_AFP=[GLOB_STRUCT.data.multiDir.ff_first_PBS] - [GLOB_STRUCT.data.multiDir.ff_first_MUSC];
bidir_Second_AFP=[GLOB_STRUCT.data.multiDir.ff_second_PBS] - [GLOB_STRUCT.data.multiDir.ff_second_MUSC];

lt_regress(bidir_Second_AFP, bidir_First_AFP, 1, 0, 1, 1, 'b');
lt_plot_zeroline;
lt_plot_zeroline_vert;


% ==== ask whether slopes are different
% - combine all first targ
X=[bidir_First_MUSC bidir_First_AFP]; % first targ
% - combine all second targ
Y = [bidir_Second_MUSC bidir_Second_AFP]; % scond targ
% -- group
group=[ones(1, length(bidir_First_MUSC)) zeros(1, length(bidir_First_AFP))]; % 1= MUSC, 0 = AFP.

aoctool(X, Y, group);





