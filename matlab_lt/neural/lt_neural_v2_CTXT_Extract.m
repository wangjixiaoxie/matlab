function [CLASSES, prms] = lt_neural_v2_CTXT_Extract(SummaryStruct, strtype)
%% lt 8/12/17 - extract pairs of contexts that are similar in any user-defined way

% outputs, for each neuron, every "branch" points that passes user-defined
% criteria

% outputs:
% regexp strings (motifs)
% sample size
% syllable to lock to (i.e. branch point);
% predur and postdur to get entire motif of interest

minN = 8;
CLASSES = struct;

%%
prms.Extract.strtype = strtype;
prms.minN = minN;



%% -- decipher what strings wanted

numsylsneeded = length(strfind(strtype, 'a'));
numvariablesyls = length(strfind(strtype, 'x'));

%% [ab]xy (i.e. convergent, followed by shared)

numbirds = length(SummaryStruct.birds);
regexprstr = [];

for i=1:numbirds
    
    numneurons = length(SummaryStruct.birds(i).neurons);
    birdname = SummaryStruct.birds(i).birdname;
    CLASSES.birds(i).birdname = birdname;
    
    for ii=1:numneurons
        
        % ========= extract song dat
        if isfield(SummaryStruct.birds(i).neurons(ii), 'dirname')
            % then is my data
        cd(SummaryStruct.birds(i).neurons(ii).dirname);
        cd ..
        
        batchf=SummaryStruct.birds(i).neurons(ii).batchfilename;
        channel_board=SummaryStruct.birds(i).neurons(ii).channel;
        [SongDat, ~, Params] = lt_neural_ExtractDat(batchf, channel_board, 0);
        else
            % then is sam/mel RA data
        [SongDat, ~, Params] = lt_neural_RASamMel_ExtractDat(SummaryStruct, i, ii);            
        end
        
        % ------------
        singlesyls = unique(SongDat.AllLabels);
        singlesyls = singlesyls(~(singlesyls=='-')); % remove dashes
        
        bb = 0;
        
        % ===================== 1 syl needed
        if numsylsneeded==1
            for j=1:length(singlesyls)
                syl1 = singlesyls(j);
                    
                    % === which string type?
                    if strtype == 'xax';
                        wildcardpositions = [1 3];
                        regexprstr = ['[a-z]' syl1 '[a-z]'];
                    end
                    
                    %%
                    
                    [start, match] = regexp(SongDat.AllLabels, regexprstr, 'start', 'match');
                    
                    [match, inds] = sort(match);
                    start = start(inds); % sort, for grp stats to be in order
                    
                    
                    
                    % ===== old version
                    N =  grpstats(start', match', 'numel');
                    matchclasses = unique(match);
                    
                    % --- remove classes with low N
                    toremove = N<minN;
                    N(toremove) = [];
                    matchclasses(toremove) = [];
                    
                    
                    
                    
                    if (0)
                    % =========== at positions with '[a-z]' make sure
                    % different across all classes
                    matchtab = tabulate(match);
                    [~, indstmp] = sort(cell2mat(matchtab(:,2)), 'descend');
                    matchtab = matchtab(indstmp,:);

                    % only keep those with >N samps;
                     matchtab = matchtab(cell2mat(matchtab(:,2)) >= minN, :);
                     
                   % --- are there >1 class?
                    if size(matchtab,1) <2
                        continue
                    end
                    
                    % ----- specific to cases with more than one wildcard
                    % STUCK HERE!!!!!!
                    wildcardsyls = nan(size(matchtab,1),length(wildcardpositions));
                    for nn=1:size(matchtab,1)
                        wildcardsyls(nn,1:length(wildcardpositions)) ...
                            = matchtab{nn,1}(wildcardpositions)
                    end
                    
                    sylstokeep = [1]; % start by keeping first one
                    if numvariablesyls >1
                        
                        % go through 
                        for nn = 2:size(matchtab,1)
                            % start from most frequent
                            wildcardsyls(sylstokeep,1)
                        end
                    end
                    end
                    
                    % =========== SAVE OUTPUT
                    disp([birdname '-neur' num2str(ii) ' found: '])
                    disp(matchclasses)
                    disp(N')
                    
                    bb = bb+1;
                    CLASSES.birds(i).neurons(ii).branchnum(bb).regexprstr = regexprstr;
                    CLASSES.birds(i).neurons(ii).branchnum(bb).matchclasses = matchclasses;
                    CLASSES.birds(i).neurons(ii).branchnum(bb).matchclasses_N = N;
            end
        end
        
        % ==================== 2 syls needed
        if numsylsneeded==2
            for j=1:length(singlesyls)
                syl1 = singlesyls(j);
                
                for jj=1:length(singlesyls)
                    syl2 = singlesyls(jj);
                    
                    % === which string type?
                    if strtype == 'xaa';
                        regexprstr = ['[a-z]' syl1 syl2];
                    elseif strtype == 'aax';
                        regexprstr = [syl1 syl2 '[a-z]'];
                    elseif strtype =='axa';
                        regexprstr = [syl1 '[a-z]' syl2];
                    end
                    
                    %%
                    
                    [start, match] = regexp(SongDat.AllLabels, regexprstr, 'start', 'match');
                    
                    [match, inds] = sort(match);
                    start = start(inds); % sort, for grp stats to be in order
                    
                    N =  grpstats(start', match', 'numel');
                    matchclasses = unique(match);
                    
                    % --- remove classes with low N
                    toremove = N<minN;
                    N(toremove) = [];
                    matchclasses(toremove) = [];
                    
                    % --- are there >1 class?
                    if length(N) <2
                        continue
                    end
                    
                    % =========== SAVE OUTPUT
                    disp([birdname '-neur' num2str(ii) ' found: '])
                    disp(matchclasses)
                    disp(N')
                    
                    bb = bb+1;
                    CLASSES.birds(i).neurons(ii).branchnum(bb).regexprstr = regexprstr;
                    CLASSES.birds(i).neurons(ii).branchnum(bb).matchclasses = matchclasses;
                    CLASSES.birds(i).neurons(ii).branchnum(bb).matchclasses_N = N;
                end
            end
        end
        % ------------------------------------
        
        % ==================== 3 syls needed
        if numsylsneeded==3
            for j=1:length(singlesyls)
                syl1 = singlesyls(j);
                
                for jj=1:length(singlesyls)
                    syl2 = singlesyls(jj);
                    
                    for jjj=1:length(singlesyls)
                        syl3 = singlesyls(jjj);
                        
                        % === which string type?
                        if strtype == 'aaxa';
                            regexprstr = [syl1 syl2 '[a-z]' syl3];
                        elseif strtype == 'xaaa'
                            regexprstr = ['[a-z]' syl1 syl2 syl3];
                        elseif strtype == 'aaax'
                            regexprstr = [syl1 syl2 syl3 '[a-z]' ];
                        end
                        
                        %%
                        
                        [start, match] = regexp(SongDat.AllLabels, regexprstr, 'start', 'match');
                        
                        [match, inds] = sort(match);
                        start = start(inds); % sort, for grp stats to be in order
                        
                        N =  grpstats(start', match', 'numel');
                        matchclasses = unique(match);
                        
                        % --- remove classes with low N
                        toremove = N<minN;
                        N(toremove) = [];
                        matchclasses(toremove) = [];
                        
                        % --- are there >1 class?
                        if length(N) <2
                            continue
                        end
                        
                        % =========== SAVE OUTPUT
                        disp([birdname '-neur' num2str(ii) ' found: '])
                        disp(matchclasses)
                        disp(N')
                        
                        bb = bb+1;
                        CLASSES.birds(i).neurons(ii).branchnum(bb).regexprstr = regexprstr;
                        CLASSES.birds(i).neurons(ii).branchnum(bb).matchclasses = matchclasses;
                        CLASSES.birds(i).neurons(ii).branchnum(bb).matchclasses_N = N;
                        
                    end
                end
            end
        end
        % ------------------------------------
        
        
        
    end
end



%% === plot mean fr differences for all "branch points"

