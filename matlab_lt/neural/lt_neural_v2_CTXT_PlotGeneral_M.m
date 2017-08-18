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
                    maxdur = params.prms.motifpredur + params.prms.ClassGeneral.frtimewindow(2);
                    numtrials = 5;
                    SylContours = [];
                    Ctxtnum = [];
                    for l = 1:numctxts
                        [sylcon, xtimes] = lt_neural_v2_ANALY_GetSylContours(dattmp.SEGEXTRACT.classnum(l).SegmentsExtract, ...
                            numtrials, maxdur);
                        
                        SylContours = [SylContours; sylcon];
                        Ctxtnum = [Ctxtnum; l*ones(size(sylcon,1),1)];
                        
                    end
                    
                    %                   AllBranchSylContours = [AllBranchSylContours SylContours];
                    %                   AllBranchSylContours_ctxtnums = [AllBranchSylContours_ctxtnums Ctxtnum];
                    
                    % ===================== COLLECT
                    bb=bb+1;
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).AllNumCtxts = numctxts;
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).AllAccuracy = accuracy;
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).AllSensitivity = sensitivity;
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).AllPiYActualMean = piYactualmean;
                    
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).AllNeurnum = kk;
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).AllBirdnum = k;
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).AllBranchnum = kkk;
                    
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).AllBranchCtxts = ctxts;
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).AllBranchRegexp = regexprstr;
                    
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).AllBranchSylContours = SylContours;
                    CLASSEScompiled.analynum(i).iterationnum(j).allbranches(bb).AllBranchSylContours_ctxtnums = Ctxtnum;
                    
                    
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



