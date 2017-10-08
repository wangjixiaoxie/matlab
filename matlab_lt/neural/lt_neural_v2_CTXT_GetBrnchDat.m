function CLASSES = lt_neural_v2_CTXT_GetBrnchDat(CLASSES, SummaryStruct, prms)
%% lt 8/12/17 - takes output of lt_neural_v2_CTXT_Extract, and gets raw dat

FFparams.collectFF=1; % note, will try to collect FF for each motif inputed in the cell array. will
FFparams.FF_PosRelToken=0; % syl to get FF of, relative to token (i.e. -1 is 1 before token;
% +1 is 1 after token
FFparams.FF_sylName=''; % Optional: what syl do you expect this to be? if incompatible will raise error
collectWNhit = 0;

% --- defaults -- don't change
                LearnKeepOnlyBase = 1;

%%

numbirds = length(CLASSES.birds);

for i=1:numbirds
   
    numneurons = length(CLASSES.birds(i).neurons);
    
    for ii=1:numneurons
        disp(['bird' num2str(i) ', neur' num2str(ii)]);
        numbranches = length(CLASSES.birds(i).neurons(ii).branchnum);
        
        if isfield(SummaryStruct.birds(i).neurons(ii), 'isRAsobermel');
            % sam/mel data
        [SongDat, NeurDat, Params] = lt_neural_RASamMel_ExtractDat(SummaryStruct, i, ii);
        LearnKeepOnlyBase=0;
        else
            % my data
        [SongDat, NeurDat, Params] = lt_neural_ExtractDat2(SummaryStruct, i, ii);
        end
                
        for iii=1:numbranches
           
            matchclasses = CLASSES.birds(i).neurons(ii).branchnum(iii).matchclasses;
            
            for mm = 1:length(matchclasses)
               
                mclass = matchclasses{mm};
                
                % ---- format to extract regexp
                n = prms.alignWhichSyl;
                mclass = [mclass(1:n-1) '(' mclass(n) ')' mclass(n+1:end)];
                
                
                % --- extract
                [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
                mclass, prms.motifpredur, prms.motifpostdur, prms.alignOnset, '', FFparams, ...
                0, 1, collectWNhit, 0, LearnKeepOnlyBase, prms.preAndPostDurRelSameTimept);


                % ================== GET FR
                if (0)
                % takes up too much space...
                    clustnum = SummaryStruct.birds(i).neurons(ii).clustnum;
                    SegmentsExtract = lt_neural_SmoothFR(SegmentsExtract, clustnum);
                end
                
                % ---- STORE
                CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT.classnum(mm).regexpstr = mclass;
                CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT.classnum(mm).SegmentsExtract = SegmentsExtract;
                CLASSES.birds(i).neurons(ii).branchnum(iii).SEGEXTRACT.classnum(mm).Params = Params;
                
            end
        end
    end
end




