function CLASSES = lt_neural_v2_CTXT_GetBrnchPosControl(CLASSES, SummaryStruct, prms, strtype)
%% lt 8/16/17 - extracts regexp for positive controls - finds 1 (rendom) set for each actual branch
% for a given strtype, code in the matching "positive control"
% CURRENT IMPLEMENTATION - gets all possible trigrams (e.g. if strtype is
% length 3). Does not yet make sure that syl is different ...

%  CURRENT STATUS - see current implementation, this is good so far
% TO DO: make sure no overlapping syls (i.e. at align syl location)
% [DONE!!]


%% === figure out what controls to use

regexpcontrol = [];

switch length(strtype)
    case 3
        regexpcontrol = '[a-z](?=[a-z][a-z])';
    case 4
        regexpcontrol = '[a-z](?=[a-z][a-z][a-z])';
end

%% extract params
% --- defaults -- don't change
LearnKeepOnlyBase = 1;
FFparams.collectFF=1; % note, will try to collect FF for each motif inputed in the cell array. will
FFparams.FF_PosRelToken=0; % syl to get FF of, relative to token (i.e. -1 is 1 before token;
% +1 is 1 after token
FFparams.FF_sylName=''; % Optional: what syl do you expect this to be? if incompatible will raise error
collectWNhit = 0;

%%

numbirds = length(CLASSES.birds);

for i=1:numbirds
    numneurons = length(CLASSES.birds(i).neurons);
    birdname = CLASSES.birds(i).birdname;
    for ii=1:numneurons
        
        disp([num2str(i) '-' num2str(ii)]);
        numbranches = length(CLASSES.birds(i).neurons(ii).branchnum);
        
       
        % --- prepare list of syls for this neuron
        [SongDat, NeurDat, Params] = lt_neural_ExtractDat2(SummaryStruct, i, ii);
        
        % ----- only keep songs without WN, if is learning expt
        if LearnKeepOnlyBase==1 & ~isfield(SummaryStruct.birds(i).neurons(ii), 'isRAsobermel'); 
            exptname = SummaryStruct.birds(i).neurons(ii).exptID;
            [islearning, LearnSummary, switchtime] = ...
                lt_neural_v2_QUICK_islearning(birdname, exptname, 1);
            
            if islearning==1
                % --- pare down labels to just those before learning start.
                songstokeep = find([NeurDat.metaDat.song_datenum]'<switchtime);
                
                if isempty(songstokeep)
                    % then no data for this neuron...
                    for iii=1:length(CLASSES.birds(i).neurons(ii).branchnum)
                        assert(~isfield(CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT.classnum(1).SegmentsExtract, 'FF_val'), 'asfasd');
                    end
                    continue
                end
                rendstokeep = find(SongDat.AllSongNum <= max(songstokeep));
                SongDat = lt_structure_subsample_all_fields(SongDat, rendstokeep);
            end
        end
        
        sylsunique = unique(SongDat.AllLabels);
        sylsunique(strfind(sylsunique, '-'))  =[];
        
        for iii=1:numbranches
            
            if ~isfield(CLASSES.birds(i).neurons(ii).branchnum(iii), 'SEGEXTRACT')
                continue
            end
            
            % -------------- 1) extract actual classes
            numclasses = length(CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT.classnum);
            motifs = {};
            Nall = [];
            
            for j=1:numclasses
                %                 if ~isfield(CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT.classnum(j).SegmentsExtract, 'FF_val')
                %                     continue
                %                 end
                
                if ~isfield(CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT.classnum(j).SegmentsExtract, 'fs')
                    continue
                end
                
                mtif = CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT.classnum(j).regexpstr;
                n = length(CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT.classnum(j).SegmentsExtract);
                
                motifs = [motifs mtif];
                Nall = [Nall n];
            end
            
            % ------------------ 1b - exit if no data...
            if isempty(Nall)
                continue
            end
            
            % ------------------ sort from large to small (so select data
            % for large first, as fewer cases can satisfy)
            [~, indstmp] = sort(Nall, 'descend');
            Nall = Nall(indstmp);
            motifs = motifs(indstmp);
            
            
            % ==================== GET MATCHING SAMPLE THAT IS POSITVE
            % CONTROL ===================================
            
            % ================ method 1 - get controls with matched
            % syllables
            
            
            % ================= method 2 - get controls, not necessarily
            % matched
            
            % ------------------------- 1) extract all potential motifs
            [start] = regexp(SongDat.AllLabels, regexpcontrol, 'start');
            
            
            % - get match syls
            strlength = length(strtype);
            indmat = [];
            for j=1:strlength
                indmat = [indmat start'+j-1];
            end
            allmotifs = tabulate(SongDat.AllLabels(indmat));
            
            
            % --------------------- 2) for each actual motif, figure out match
            repickmotifs =1; % will keep going until get motifs that don't have overlapping sysl.
            counter=1;
            while repickmotifs==1
                allmotifs_tmp = allmotifs; % because will pare down ...
                MotifsPicked = {};
                
                for j=1:length(motifs)
                    
                    n = Nall(j); % need this sample size
                    inds_potential = find(cell2mat(allmotifs_tmp(:,2))>=n); % those with enough sample size
                    assert(~isempty(inds_potential), 'asdfsadf huh?');
                    indtmp = inds_potential(randi(length(inds_potential), 1)); % pick a random one
                    
                    MotifsPicked{j} = allmotifs_tmp{indtmp, 1};
                    
                    % -- remove this motif from consideration
                    allmotifs_tmp(indtmp,:) = [];
                end
                
                % make sure syls (for each ordinal position) are different
                % between all classes.
                tmp = cell2mat(MotifsPicked');
                repickmotifs = 0;
                for j=1:size(tmp,2)
                    if length(unique(tmp(:,j))) ~= length(tmp(:, j))
                        % then some syls are the samne ..
                        repickmotifs = 1;
                    end
                end
                counter = counter+1;
                if counter==100
                    disp('NOTE: failed to get totally nonoverlapping motifs ...');
                    break
                end
            end
            
            % ------------------- 3) extract regexp for these motifs
            for j=1:length(MotifsPicked)
                mclass = MotifsPicked{j};
                
                % ---- put parantheses for token
                mclass = [mclass(1:prms.alignWhichSyl-1) '(' ...
                    mclass(prms.alignWhichSyl) ')' mclass(prms.alignWhichSyl+1:end)];
                
                % --- extract
                if  ~isfield(SummaryStruct.birds(i).neurons(ii), 'isRAsobermel')
                    learntmp = LearnKeepOnlyBase;
                else
                    learntmp = 0;
                end
                [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
                    mclass, prms.motifpredur, prms.motifpostdur, prms.alignOnset, '', FFparams, ...
                    0, 1, collectWNhit, 0, learntmp, prms.preAndPostDurRelSameTimept);
                
                % --- confirm that sample sizes match actual data
                assert(length(SegmentsExtract) >= Nall(j), 'problem - likely becuyase regexp didnt get overlapping, so N is small')
                
                % --- downsample to size of actual data(remove random dat)
                SegmentsExtract(randperm(length(SegmentsExtract), length(SegmentsExtract)-Nall(j))) =[];
                
                
                % ---- GET FR
                if (0)
                clustnum = SummaryStruct.birds(i).neurons(ii).clustnum;
                SegmentsExtract = lt_neural_SmoothFR(SegmentsExtract, clustnum);
                end
                
                
                % ===================== SAVE IN OUTPUT
                CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT_POSCONTR.classnum(j).regexpstr = mclass;
                CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT_POSCONTR.classnum(j).SegmentsExtract = SegmentsExtract;
                CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT_POSCONTR.classnum(j).Params = Params;
            end
        end
    end
end

