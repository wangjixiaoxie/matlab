function ALLBRANCH = lt_neural_v2_CTXT_BranchGetFR(ALLBRANCH, saveOn)

%% lt 8/27/17 - get smoothed FR for all branches... and saves

FFparams.collectFF=1; % note, will try to collect FF for each motif inputed in the cell array. will
FFparams.FF_PosRelToken=0; % syl to get FF of, relative to token (i.e. -1 is 1 before token;
% +1 is 1 after token
FFparams.FF_sylName=''; % Optional: what syl do you expect this to be? if incompatible will raise error
collectWNhit = 0;
LearnKeepOnlyBase = 1;

%%  get max neurons, and branches

Maxneurons = [];
Maxbranches = [];

numalign = length(ALLBRANCH.alignpos);

for i=1:numalign
    
    numbirds = length(ALLBRANCH.alignpos(i).bird);
    
    for ii = 1:numbirds
        numbranches = length(ALLBRANCH.alignpos(i).bird(ii).branch);
        
        Maxbranches = max([Maxbranches numbranches]);
        
        for bb = 1:numbranches
            
            numneurons = length(ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron);
            
            Maxneurons = max([Maxneurons numneurons]);
            
        end
        
    end
end

%% ==================================== PLOT SMOOTHED FR DPRIME (AS IN SOBER/WOHLGEMUTH)
numalign = length(ALLBRANCH.alignpos);
DprimeAll = [];

for i=1:numalign
%     alignsyl = ALLBRANCH.alignpos(i).alignsyl;
%     alignons = ALLBRANCH.alignpos(i).alignonset;
    
    numbirds = length(ALLBRANCH.alignpos(i).bird);
    
    %     Xall = [];
    %     Yall = [];
    %     Xcont = [];
    %     Ycont = [];
    %
    %     Xcell = {};
    %     Ycell = {};
    %     Xcontcell = {};
    %     Ycontcell = {};
    
    for ii=1:numbirds
        
        for nn = 1:Maxneurons
            
            % ------- if data exist, then extract song dat
            if length(ALLBRANCH.SummaryStruct.birds(ii).neurons) <nn
                continue
            end
            
            % ------ extract dat
            [SongDat, NeurDat, Params] = lt_neural_ExtractDat2(ALLBRANCH.SummaryStruct, ii, nn);
            
            for bb = 1:Maxbranches
                if length(ALLBRANCH.alignpos(i).bird(ii).branch)<bb
                    continue
                end
                if length(ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron)<nn
                    continue
                end
                if isempty(ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).yvals)
                    continue
                end
                
                
                % =============================== ACTUAL DAT
                regexpstrlist = ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).prms_regexpstrlist;
                
                numclasses = length(regexpstrlist);
                assert(numclasses>1,'sdfasd');
                
                motifpredur = ALLBRANCH.alignpos(i).ParamsFirstIter.motifpredur;
                motifpostdur = ALLBRANCH.alignpos(i).ParamsFirstIter.motifpostdur;
                alignonset = ALLBRANCH.alignpos(i).ParamsFirstIter.alignOnset;
                preAndPostDurRelSameTimept = ALLBRANCH.alignpos(i).ParamsFirstIter.preAndPostDurRelSameTimept;
                assert(alignonset == ALLBRANCH.alignpos(i).alignonset);
                
                Ndat = [];
                for cc = 1:numclasses
                    
                    regexpstr = regexpstrlist{cc};
                    
                    % ----- for this branch, collect segextract
                    [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
                        regexpstr, motifpredur, motifpostdur, alignonset, '', FFparams, ...
                        0, 1, collectWNhit, 0, LearnKeepOnlyBase, preAndPostDurRelSameTimept);
                    
                   Ndat = [Ndat length(SegmentsExtract)];
                    
                    % ----- collect smoothed FR
                    clustnum = ALLBRANCH.SummaryStruct.birds(ii).neurons(nn).clustnum;
                    SegmentsExtract = lt_neural_SmoothFR(SegmentsExtract, clustnum);
                    
                    % ============= save smoothed FR
                    ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).FR.classnum(cc).regexpstr = regexpstr;
                    ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).FR.classnum(cc).FRsmooth_rate_CommonTrialDur ...
                        = single([SegmentsExtract.FRsmooth_rate_CommonTrialDur]);
                    ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).FR.classnum(cc).FRsmooth_xbin_CommonTrialDur ...
                        = single(SegmentsExtract(1).FRsmooth_xbin_CommonTrialDur);
                end
                
                
                % ================================== POS CONTROL IN
                % PROGRESS ~~~~
                % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                regexpstrlist = ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).prms_regexpstrlist_POSCONTR;
                numclasses = length(regexpstrlist);
                Nall = ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).prms_POSCONTR_regexpstrlist_N;
                
                Ndat_pos = [];
                for cc = 1:numclasses
                                        
                    regexpstr = regexpstrlist{cc};
                    
                    % ----- for this branch, collect segextract
                    [SegmentsExtract, Params]=lt_neural_RegExp(SongDat, NeurDat, Params, ...
                        regexpstr, motifpredur, motifpostdur, alignonset, '', FFparams, ...
                        0, 1, collectWNhit, 0, LearnKeepOnlyBase, preAndPostDurRelSameTimept);

                    
                     % --- downsample to size of actual data(remove random dat)
                SegmentsExtract(randperm(length(SegmentsExtract), length(SegmentsExtract)-Nall(cc))) =[];

                Ndat_pos = [Ndat_pos length(SegmentsExtract)];
                                   
                    % ----- collect smoothed FR
                    clustnum = ALLBRANCH.SummaryStruct.birds(ii).neurons(nn).clustnum;
                    SegmentsExtract = lt_neural_SmoothFR(SegmentsExtract, clustnum);
                    
                    % =============== save
                    ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).FR_POSCONTR.classnum(cc).regexpstr = regexpstr;
                    ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).FR_POSCONTR.classnum(cc).FRsmooth_rate_CommonTrialDur ...
                        = single([SegmentsExtract.FRsmooth_rate_CommonTrialDur]);
                    ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).FR_POSCONTR.classnum(cc).FRsmooth_xbin_CommonTrialDur ...
                        = single([SegmentsExtract(1).FRsmooth_xbin_CommonTrialDur]);

                end
               
                if (0)
                    % debugging
                disp([num2str(Ndat) ' -- ' num2str(Ndat_pos)]);
               
                pause;
                
                end
                
            end
        end
    end
end
    
%% ==============
if saveOn ==1
    % --- go to save dir
    savedir = '/bluejay5/lucas/analyses/neural/CTXT_ClassGeneral_M';
    cd(savedir)
    
    % --- save
    tstamp = lt_get_timestamp(0);
    strtype = ALLBRANCH.alignpos(1).ParamsFirstIter.Extract.strtype;
    savefname = ['ALLBRANCH_' strtype '_' tstamp];
    
    save(savefname, 'ALLBRANCH');
    
end

