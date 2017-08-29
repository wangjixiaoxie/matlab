function [CLASSEScompiled] = lt_neural_v2_CTXT_PlotGeneral_M(strtype, algnsyl, algnonset)
%% lt 8/15/17 -

% strtype = 'xaaa' - i.e. what was used for lt_neural_v2_CTXT_Extract. Will
% go thru all saved analysis folders. any that match will be plotted

% SummaryStruct - should be same as the one used for this analysis. Will
% loook in folder for saved version. will only use this arguemnt if that
% saved ver. doesn't exist.


savedir = '/bluejay5/lucas/analyses/neural/CTXT_ClassGeneral_M';

%%

cd(savedir)

%%

CLASSEScompiled = struct;

dirnames = dir(['Results_' strtype '_AlgnSyl' num2str(algnsyl) 'Onset' num2str(algnonset) '_*']);

for i=1:length(dirnames)
    
    uscores = strfind(dirnames(i).name, '_');
    
    analydate = dirnames(i).name(uscores(3)+1:uscores(5)-1);
    savenote = dirnames(i).name(uscores(5)+1:end);
    
    cd(dirnames(i).name);
    disp('-------------------');
    
    % =========== LOAD ALL SAVED STRUCTURES, AND COLLECT DATA ACROSS ALL
    datfiles = dir('classes*');
    
    % -- sort by date
    [~, inds] = sort([datfiles.datenum]);
    datfiles = datfiles(inds);
    
    % ==== collect stuff
    CLASSEScompiled.analynum(i).analydate = analydate;
    CLASSEScompiled.analynum(i).savenote = savenote;
    
    tmp = load('SummaryStruct');
    CLASSEScompiled.analynum(i).SummaryStruct = tmp.SummaryStruct;
    
    for j=1:length(datfiles)
        
        dotind = strfind(datfiles(j).name, '.mat');
        
        filenum = datfiles(j).name(8:dotind(1)-1);
        
        %         disp(filenum)
        classes = load(['classes' num2str(filenum) '.mat']);
        params = load(['params' num2str(filenum) '.mat']);
        
        
        % ==== collect information about this analysis instance
        CLASSEScompiled.analynum(i).iterationnum(j).filenum = filenum;
        CLASSEScompiled.analynum(i).iterationnum(j).params = params.prms;
        
        
        %         AllNumCtxts = [];
        %         AllAccuracy = [];
        %         AllSensitivity = [];
        %
        %         AllNeurnum = [];
        %         AllBirdnum = [];
        %         AllBranchnum = [];
        %
        %         AllBranchCtxts = {};
        %         AllBranchRegexp = {};
        %
        %         AllBranchSylContours = {};
        %         AllBranchSylContours_ctxtnums = {};
        
        bb=0; % branch counter;
        disp([num2str(i) '-' num2str(j)]);
        % ==== collect dat
        numbirds = length(classes.CLASSES.birds);
        for k=1:numbirds
            numneurons = length(classes.CLASSES.birds(k).neurons);
            
            for kk=1:numneurons
                numbranches = length(classes.CLASSES.birds(k).neurons(kk).branchnum);
                
                for kkk=1:numbranches
                    
                    dattmp = classes.CLASSES.birds(k).neurons(kk).branchnum(kkk);
                    
                    if ~isfield(dattmp, 'CLASSIFIER')
                        continue
                    end
                    
                    if isempty(dattmp.CLASSIFIER)
                        continue
                    end
                    
                    if isempty(dattmp.CLASSIFIER.Rslts_accuracy)
                        continue
                    end
                    
                    % --
                    numctxts = length(dattmp.CLASSIFIER.ctxts_that_exist);
                    accuracy = dattmp.CLASSIFIER.Rslts_accuracy;
                    sensitivity = dattmp.CLASSIFIER.Rslts_sensitivity_mean;
                    piYactualmean = [];
                    if isfield(dattmp.CLASSIFIER, 'Rslts_PiYActual')
                        piYactualmean = mean(dattmp.CLASSIFIER.Rslts_PiYActual);
                    end
                       
                    
                    ctxts = dattmp.CLASSIFIER.ctxts_that_exist_name;
                    regexprstr = dattmp.regexprstr;
                    
                    % --- out
                    %                   AllNumCtxts = [AllNumCtxts numctxts];
                    %                   AllAccuracy = [AllAccuracy accuracy];
                    %                   AllSensitivity = [AllSensitivity sensitivity];
                    %
                    %                   AllNeurnum = [AllNeurnum kk];
                    %                   AllBirdnum = [AllBirdnum k];
                    %                   AllBranchnum = [AllBranchnum kkk];
                    %
                    %                   AllBranchCtxts = [AllBranchCtxts ctxts];
                    %                   AllBranchRegexp = [AllBranchRegexp regexprstr];
                    
                    
                    % ---- extract syl contours
%                     maxdur = params.prms.motifpredur + params.prms.ClassGeneral.frtimewindow(2);
                    maxdur = params.prms.motifpredur + params.prms.motifpostdur;
                    
                    numtrials = 4;
                    SylContours = [];
                    Ctxtnum = [];
                    numctxtstmp = length(dattmp.SEGEXTRACT.classnum);
                    for l = 1:numctxtstmp
                        
                        if ~isfield(dattmp.SEGEXTRACT.classnum(l).SegmentsExtract, 'fs');
                            continue
                        end
                        
                        [sylcon, xtimes] = lt_neural_v2_ANALY_GetSylContours(dattmp.SEGEXTRACT.classnum(l).SegmentsExtract, ...
                            numtrials, maxdur);
                        
                        SylContours = [SylContours; sylcon];
                        Ctxtnum = [Ctxtnum; l*ones(size(sylcon,1),1)];
                    end
                    assert(sum(sum(isnan(SylContours)))./numel(SylContours) < 0.01, 'problem - >1% are nan');
                    SylContours = int8(SylContours);
                    %                   AllBranchSylContours = [AllBranchSylContours SylContours];
                    %                   AllBranchSylContours_ctxtnums = [AllBranchSylContours_ctxtnums Ctxtnum];
                    
                    % -- CONFMAT
                    ConfMat = dattmp.CLASSIFIER.Rslts_ConfMat;
                    
                    
                    % ========== CONTROLS
                    if isfield(dattmp.CLASSIFIER, 'ShuffNeg_ConfMatAll')
                    NEGCONTR_ConfMatAll = dattmp.CLASSIFIER.ShuffNeg_ConfMatAll;
                    POSTCONTR_ConfMat = dattmp.CLASSIFIER.ContrPos_ConfMat;
                    end

                    
                    % ===================== COLLECT
                    bb=bb+1;
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).AllNumCtxts = numctxts;
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).AllAccuracy = accuracy;
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).AllSensitivity = sensitivity;
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).AllPiYActualMean = piYactualmean;
                    
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).AllNeurnum = kk;
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).AllBirdnum = k;
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).AllBranchnum = kkk;
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).AllBranchnum_NOTE = 'AllBranchNum is meaningless (is arbitrary, originally was lower level than each neuron, used in CLASSES). Also, the index for allbranches is also meaningless';
                    
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).AllBranchCtxts = ctxts;
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).AllBranchRegexp = regexprstr;
                    
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).AllBranchSylContours = SylContours;
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).AllBranchSylContours_ctxtnums = Ctxtnum;
                    
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).AllConfMat = int16(ConfMat);
                    
                    % ========== COLLECT FR BINNED THAT WAS USED FOR
                    % CLASSIFICATION
                    % -- dat
                    frtbins = dattmp.CLASSIFIER.Dat_xbins;
                    frmean = mean(dattmp.CLASSIFIER.Dat_X,1);
                    frstd_xtrials = std(dattmp.CLASSIFIER.Dat_X,0,1);
                    frstd_xbins = median(std(dattmp.CLASSIFIER.Dat_X, 0, 2));
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).FRdat_t = single(frtbins);
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).FRdat_mean = single(frmean);
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).FRdat_stdxtrials = single(frstd_xtrials);
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).FRdat_stdxbins = single(frstd_xbins);
                    
                    % -- pos contr
                    frmean = mean(dattmp.CLASSIFIER.ContrPos_Xall,1);
                    frstd_xtrials = std(dattmp.CLASSIFIER.ContrPos_Xall,0,1);
                    frstd_xbins = median(std(dattmp.CLASSIFIER.ContrPos_Xall,0,2));
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).FRdat_mean_PosContr = single(frmean);
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).FRdat_stdxtrials_PosContr = single(frstd_xtrials);
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).FRdat_stdxbins_PosContr = single(frstd_xbins);
                   
                    
                    % --- CONTROLS
                    if isfield(dattmp.CLASSIFIER, 'ShuffNeg_ConfMatAll')
                    % if has multiple, then take random one (and save
                    % multiple as backup)
                    if length(size(NEGCONTR_ConfMatAll)) >2
                        % then take random one
                        negconfmat = NEGCONTR_ConfMatAll(:,:, randi(size(NEGCONTR_ConfMatAll,3)));
                    else
                        negconfmat = NEGCONTR_ConfMatAll;
                    end
                    
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).NEGCONTR_AllConfMatAllShuff = int16(NEGCONTR_ConfMatAll);
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).NEGCONTR_AllConfMat = int16(negconfmat);
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).POSCONTR_AllConfMat = int16(POSTCONTR_ConfMat);
                    end                   
                    
                    % --------------- save segextract
                    if (0) % TOO LARGE
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).SEGEXTRACT = ...
                        dattmp.SEGEXTRACT;

                    numclasses= length(dattmp.SEGEXTRACT.classnum);
                    for cc=1:numclasses
                        numtrials = length(CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).SEGEXTRACT.classnum(cc).SegmentsExtract);
                        
                        for ccc = 1:numtrials
                        CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).SEGEXTRACT.classnum(cc).SegmentsExtract(ccc).spk_Clust = ...
                            int8(CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).SEGEXTRACT.classnum(cc).SegmentsExtract(ccc).spk_Clust);
                        CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).SEGEXTRACT.classnum(cc).SegmentsExtract(ccc).spk_Times = ...
                            single(CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).SEGEXTRACT.classnum(cc).SegmentsExtract(ccc).spk_Times);
                        CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).SEGEXTRACT.classnum(cc).SegmentsExtract(ccc).sylOnTimes_RelDataOnset = ...
                            single(CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).SEGEXTRACT.classnum(cc).SegmentsExtract(ccc).sylOnTimes_RelDataOnset);
                        CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).SEGEXTRACT.classnum(cc).SegmentsExtract(ccc).sylOffTimes_RelDataOnset = ...
                            single(CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).SEGEXTRACT.classnum(cc).SegmentsExtract(ccc).sylOffTimes_RelDataOnset);
                        end
                    end
                    
                    % ----------------- save segextract POS CONTROL
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).SEGEXTRACT_POSCONTR = ...
                        dattmp.SEGEXTRACT_POSCONTR;

                    numclasses= length(dattmp.SEGEXTRACT_POSCONTR.classnum);
                    for cc=1:numclasses
                        numtrials = length(CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).SEGEXTRACT_POSCONTR.classnum(cc).SegmentsExtract);
                        
                        for ccc = 1:numtrials
                        CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).SEGEXTRACT_POSCONTR.classnum(cc).SegmentsExtract(ccc).spk_Clust = ...
                            int8(CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).SEGEXTRACT_POSCONTR.classnum(cc).SegmentsExtract(ccc).spk_Clust);
                        CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).SEGEXTRACT_POSCONTR.classnum(cc).SegmentsExtract(ccc).spk_Times = ...
                            single(CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).SEGEXTRACT_POSCONTR.classnum(cc).SegmentsExtract(ccc).spk_Times);
                        CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).SEGEXTRACT_POSCONTR.classnum(cc).SegmentsExtract(ccc).sylOnTimes_RelDataOnset = ...
                            single(CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).SEGEXTRACT_POSCONTR.classnum(cc).SegmentsExtract(ccc).sylOnTimes_RelDataOnset);
                        CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).SEGEXTRACT_POSCONTR.classnum(cc).SegmentsExtract(ccc).sylOffTimes_RelDataOnset = ...
                            single(CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).SEGEXTRACT_POSCONTR.classnum(cc).SegmentsExtract(ccc).sylOffTimes_RelDataOnset);
                        end
                    end
                    end
                    
                    % ============ save regexp strings for controls
                    POSCONTR_regexpstrlist ={dattmp.SEGEXTRACT_POSCONTR.classnum.regexpstr};
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).POSCONTR_regexpstrlist = POSCONTR_regexpstrlist;
                    
                    POSCONTR_regexpstrlist_N = cellfun('length', {dattmp.SEGEXTRACT_POSCONTR.classnum.SegmentsExtract});
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).POSCONTR_regexpstrlist_N = POSCONTR_regexpstrlist_N;
                end
                
            end
            
        end
        
    end
    
    cd ..
end

%% ==== save otuptu

tstamp = lt_get_timestamp(0);

fname = ['CLASSEScompiled_' strtype '_AlgnSyl' num2str(algnsyl) 'Onset' num2str(algnonset) '_' tstamp];

save(fname, 'CLASSEScompiled');



