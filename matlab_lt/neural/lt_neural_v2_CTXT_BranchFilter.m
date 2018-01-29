function ALLBRANCH = lt_neural_v2_CTXT_BranchFilter(ALLBRANCH, Params)
%% takes ALLBRANCH and outputs modified ALLBRANCH
% === removes branches that don't pass criteria (makes them empty)
% == IMPROTANT: bird, neuron, branch indices will all still be correct
% (e.g. same as in SummaryStruct). will simply empty it of data.

%% PARAMS

% Params.LocationsToKeep = {};
% % Params.birdstoexclude = {'bk7', 'bu77wh13', 'or74bk35', 'wh6pk36', 'br92br54'};
% Params.birdstoexclude = {}; % performs 1st
% Params.birdstokeep = {}; % performs 2nd
% Params.RemoveRepeats = 0; % Removes if, for any class in branch, presyl is same as token (e.g. a(a)a)
% Params.durThreshOmega.syl = 0.2; % omega2 (will only keep if lower) [leave empty to ignore]
% Params.durThreshOmega.gappre= [];
% Params.durThreshOmega.gappost= [];
% Params.GapDurPreMax = 0.5; % then will throw out if median pregap dur (for any
%     % class within branch) is longer than this (sec)
% Params.RemoveHandCoded =0 ; % see below
% Params.expttokeep = {'RAlearn1'};


%%
LocationsToKeep = Params.LocationsToKeep;
birdstoexclude = Params.birdstoexclude;
RemoveRepeats = Params.RemoveRepeats;
durThreshOmega = Params.durThreshOmega;
GapDurPreMax = Params.GapDurPreMax;
RemoveHandCoded = Params.RemoveHandCoded;

if ~isfield(Params, 'birdstokeep')
    birdstokeep = {};
else
    birdstokeep = Params.birdstokeep;
end

if ~isfield(Params, 'expttokeep')
    Params.expttokeep = {}; % then assumes get all expts
end

%%

apos = 1; % assumes only one analysis done


%% === HAND CODED BAD BRANCHES (bad labeling)

HandCoded(1).birdname = 'B53O71';
HandCoded(1).branches = {'[a-z]kk', '[a-z]gg', '[a-z]ca', '[a-z]aj', '[a-z]ac'};

HandCoded(2).birdname = 'Pu55Pu22';
HandCoded(2).branches = {'[a-z]ca'};

HandCoded(3).birdname = 'W32Pi51';
HandCoded(3).branches = {'[a-z]dd'};


%% extract duration (gap, syl) anovas
% ALLBRANCH = lt_neural_v2_CTXT_BranchGaps(ALLBRANCH);


%%
numbirds = length(ALLBRANCH.alignpos(apos).bird);
%     motifpredur = ALLBRANCH.alignpos(apos).ParamsFirstIter.motifpredur;
Birdstoremove = [];
for ii=1:numbirds
    
    % #################################### REMOVE BIRDS?
    birdname = ALLBRANCH.SummaryStruct.birds(ii).birdname;
    if any(strcmp(birdstoexclude, birdname))
        disp(['REMOVING ' birdname]);
        Birdstoremove = [Birdstoremove ii];
        continue
    end
    
    if ~isempty(birdstokeep)
       % --- then only continue if this bird is wanted
       if ~any(strcmp(birdstokeep, birdname))
                   disp(['REMOVING ' birdname]);
        Birdstoremove = [Birdstoremove ii];

       continue
    end
       
    end
    
    
    numbranch = length(ALLBRANCH.alignpos(apos).bird(ii).branch);
    for j=1:numbranch
        numneuron = length(ALLBRANCH.alignpos(apos).bird(ii).branch(j).neuron);
        
        Neuronstoremove = [];
        for nn=1:numneuron
            
            datneur = ALLBRANCH.alignpos(apos).bird(ii).branch(j).neuron(nn);
            if isempty(datneur.yvals)
                continue
            end
            
            % ################## EXPT
            if ~isempty(Params.expttokeep)
            expt_this = ALLBRANCH.SummaryStruct.birds(ii).neurons(nn).exptID;
            
            if ~any(strcmp(Params.expttokeep, expt_this))
                disp(['[EXPT removed] - ' expt_this]);
                Neuronstoremove = [Neuronstoremove nn];
                
                continue
            end
            
            end
            
            % #################################### LOCATION
            if ~isempty(LocationsToKeep)
                if isfield(ALLBRANCH.SummaryStruct.birds(ii).neurons(nn),'isRAsobermel')
                    location = 'RA';
                else
                    location = ALLBRANCH.SummaryStruct.birds(ii).neurons(nn).NOTE_Location;
                end
                if ~any(strcmp(LocationsToKeep, location))
                    % then skip this neuron
                    Neuronstoremove = [Neuronstoremove nn];
                    disp(['[LOCATION] removed neur ' num2str(nn)]);
                    continue
                end
            end
            
            % ################################ REPEATS?
            if RemoveRepeats ==1
                isrepeat=0;
                for k =1:length(datneur.prms_regexpstrlist)
                    thisstr = datneur.prms_regexpstrlist{k};
                    tokensyl = thisstr(strfind(thisstr, '(')+1);
                    presyl = thisstr(strfind(thisstr, '(')-1);
                    
                    if tokensyl == presyl
                        isrepeat =1;
                    end
                end
                
                if isrepeat==1
                    % remove this neuron
                    Neuronstoremove = [Neuronstoremove nn];
                    disp(['[REPEATS] removed neur ' num2str(nn)]);
                    continue
                end
            end
            
            
            % ##################################### syl/gap duration
            % differences between classes for this branch
            if datneur.DurAnovas.syl_omega>durThreshOmega.syl | ...
                    datneur.DurAnovas.gappre_omega > durThreshOmega.gappre | ...
                    datneur.DurAnovas.gappost_omega > durThreshOmega.gappost
                
                % ------------- remove neuron
                Neuronstoremove = [Neuronstoremove nn];
                disp(['[GAP/SYL DUR] removed neur ' num2str(nn)]);
                continue
            end
            
            
            % ################################## pre gap dur too long?
            if ~isempty(GapDurPreMax)
            functmp = @(x)median(x);
            pregapdurs_med = cellfun(functmp, {datneur.SylGapDurs.classnum.Dur_gappre});
            if any(pregapdurs_med > GapDurPreMax)

                % ------------- remove neuron
                Neuronstoremove = [Neuronstoremove nn];
                disp(['[PRE GAP TOO LONG] removed neur ' num2str(nn)]);
                continue

            end
            end
            
            % ####################################### HAND CODED, REMOVE
            % CERTAIN BRANCHES (BY REGEXP) - I.E. BAD LABELING
            if RemoveHandCoded==1
            if any(strcmp({HandCoded.birdname}, birdname));
                thisbranch = datneur.prms_regexpstr;
                thisind = strcmp({HandCoded.birdname}, birdname);
                if any(strcmp(HandCoded(thisind).branches, thisbranch))
                    % then remove this neuron
                    Neuronstoremove = [Neuronstoremove nn];
                    disp(['[HAND CODED] removed neur ' num2str(nn)]);
                    continue
                end
                
            end
            end
        end
        
        % ############################# REMOVE NEURONS
        for nn = Neuronstoremove;
            fnames = fieldnames(ALLBRANCH.alignpos(apos).bird(ii).branch(j).neuron(nn));
            for ff = fnames'
                ALLBRANCH.alignpos(apos).bird(ii).branch(j).neuron(nn).(ff{1}) = [];
            end
        end
    end
    
end

%% ################################### REMOVE BIRDS
for bb = Birdstoremove
    
    % --- remove all neurons
    numbranch = length(ALLBRANCH.alignpos(apos).bird(bb).branch);
    for j=1:numbranch
        
        numneuron = length(ALLBRANCH.alignpos(apos).bird(bb).branch(j).neuron);
        for nn=1:numneuron
            
            % ======================== REMOVE!!! (MAKE ALL FIELDS EMPTY)
            fnames = fieldnames(ALLBRANCH.alignpos(apos).bird(bb).branch(j).neuron(nn));
            for ff = fnames'
                ALLBRANCH.alignpos(apos).bird(bb).branch(j).neuron(nn).(ff{1}) = [];
            end
        end
    end
end

