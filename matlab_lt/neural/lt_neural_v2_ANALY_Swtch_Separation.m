function [Datstruct, getdprime] = lt_neural_v2_ANALY_Swtch_Separation(MOTIFSTATS_Compiled, ...
    SwitchStruct, onlyUseDatOnSwitchDays)
%%  lt 9/7/17 - compare separation between contexts pre and end of

Nmin = 5;
getdprime=0;

anovapostsylonsetmax = 0.075; % to cut off anova at this after syl onset.
maxduranova = MOTIFSTATS_Compiled.birds(1).exptnum(1).MOTIFSTATS.params.motif_predur ...
    + anovapostsylonsetmax;

windowanova = 0.02;
windshiftanova = 0.01;


%%
% learning
% - does not compare same neuron pre and post. instead compares pairs.
% depending on extp (1) targ one context; 2) bidir; 3) same dir, predict
% difference effect on different types of pairs.


%
%
% %% == all cases where target only one context for a given syl.
%
%
% %
% % MOTIFSTATS_Compiled.birds(i).exptnum(ii).MOTIFSTATS.params
%
% targsyls = SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).learningContingencies(1:2:end);
%
%
% %%
%

% % 1) ONSETS
% - one targ
% - two targ (bidir)
% - two targ (same dir)
%
% % 2) SWITCHES
% - any transition between (i) onetarg, (ii) two targ(bidir) and (iii) two targ (same dir)



%% ---------------------------------

Datstruct = struct;

% ---
numbirds = length(MOTIFSTATS_Compiled.birds);
ecount = 0;

for bb = 1:numbirds
    
    numexpts = length(MOTIFSTATS_Compiled.birds(bb).exptnum);
    birdname = MOTIFSTATS_Compiled.birds(bb).birdname;
    for ee =1:numexpts
        
        numswitches = length(SwitchStruct.bird(bb).exptnum(ee).switchlist);
        motifstats = MOTIFSTATS_Compiled.birds(bb).exptnum(ee).MOTIFSTATS;
        exptname = SwitchStruct.bird(bb).exptnum(ee).exptname;
        
        for ss = 1:numswitches
            
            % ===================== only continue if is onset of experiment
            tmptmp = SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).learningContingencies(2:2:end);
            tmptmp = cell2mat(tmptmp');
            if ~all(tmptmp(:,1)==0)
                disp('SKIPPED - not starting from WN off')
                continue
            end
            
            % ======= note down experiment type
            singlesyls = motifstats.params.SingleSyls_inorder;
            learnconting = SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).learningContingencies;
            
            expttype = [];
            numtargs = [];
            if length(learnconting) ==2
                numtargs = 1;
                % then only one targ - other contexts exist?
                targsinglesyl = learnconting{1}(findstr(learnconting{1}, '(')+1);
                if sum(strcmp(singlesyls, targsinglesyl))>1
                    expttype =   'one targ context';
                elseif sum(strcmp(singlesyls, targsinglesyl)) ==1
                    expttype = 'one targ - stereotyped';
                else
                    asdfhasdlhihildsaf;
                end
            elseif length(learnconting)>2
                numtargs = length(learnconting)/2;
                % then multiple targets - are they same direction?
                if length(unique(cell2mat(SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).learningDirs(2:2:end)))) ==1
                    expttype = 'mult targ context - samedir';
                else
                    expttype = 'mult targ context - diff dir';
                end
            end
            
            ecount = ecount+1;
            Datstruct.exptcount(ecount).birdnum_old = bb;
            Datstruct.exptcount(ecount).exptnum_old = ee;
            Datstruct.exptcount(ecount).switchnum_old = ss;
            Datstruct.exptcount(ecount).expttype = expttype;
            Datstruct.exptcount(ecount).numtargs = numtargs;
            
            Datstruct.exptcount(ecount).learndirs = SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).learningDirs;
            
            
            % ########################## display things about this expt
            disp(' ');
            disp(' ############################################## ');
            disp([birdname '-' exptname '-sw' num2str(ss)]);
            tmp = datestr(SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).switchdnum, 'ddmmmyyyy-HHMM');
            disp(['time of switch: ' tmp]);
            disp(['expt type: ' expttype]);
            disp(['targets: ' SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).learningDirs]);
            disp(['same type syls: ' motifstats.params.SameTypeSyls]);
            
            % #####################
            
            % ==========================================================
            % collect data for each neuron
            goodneurons = SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).goodneurons;
            
            for nn=goodneurons
                
                motiflist = motifstats.params.motif_regexpr_str;
                nummotifs = length(SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).neuron(nn).DATA.motif);
                mcount = 0;
                %                 nummotifs = length(motifstats.params.motif_regexpr_str);
                
                % ################################# display song labeling
                % progress for this neuron
                
                cd(MOTIFSTATS_Compiled.birds(bb).exptnum(ee).SummaryStruct.birds(1).neurons(nn).dirname);
                batchfile = MOTIFSTATS_Compiled.birds(bb).exptnum(ee).SummaryStruct.birds(1).neurons(nn).batchfilename;
                eval(['!cp ' batchfile ' ..']); % copy batch file one dir up. [as this is most accurate batch file]
                cd .. % labeled files are up one.
                disp(['neuron ' num2str(nn) ' out of ' num2str(max(goodneurons))]);
                disp(MOTIFSTATS_Compiled.birds(bb).exptnum(ee).SummaryStruct.birds(1).neurons(nn).dirname)
                
                lt_batch_disp_all_labels(batchfile);
                
                
                
                % ###################################
                
                for mm=1:nummotifs
                    syl1 = motiflist{mm};
                    
                    for mmm = mm+1:nummotifs
                        syl2 = motiflist{mmm};
                        
                        
                        
                        
                        % ==================================== COLLECT
                        % DPRIME
                        
                        
                        if onlyUseDatOnSwitchDays==0
                        % --------------------------------------- train and base Inds
                        baseInds_m1 = find(SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).neuron(nn).DATA.motif(mm).baseInds);
                        trainInds_m1 = find(SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).neuron(nn).DATA.motif(mm).trainInds);
                        
                        baseInds_m2 = find(SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).neuron(nn).DATA.motif(mmm).baseInds);
                        trainInds_m2 = find(SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).neuron(nn).DATA.motif(mmm).trainInds);
                        else
                        % --------------------------------------- train and base Inds
                        baseInds_m1 = find(SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).neuron(nn).DATA.motif(mm).baseInds_WithinDayOfSw);
                        trainInds_m1 = find(SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).neuron(nn).DATA.motif(mm).trainInds_WithinDayOfSw);
                        
                        baseInds_m2 = find(SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).neuron(nn).DATA.motif(mmm).baseInds_WithinDayOfSw);
                        trainInds_m2 = find(SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).neuron(nn).DATA.motif(mmm).trainInds_WithinDayOfSw);
                        end
                        
                        % --- figure out max sample size for base and train that allows for matched
                        % sample sizes
                        
                        N1 = min([length(baseInds_m1) length(trainInds_m1)]);
                        N2 = min([length(baseInds_m2) length(trainInds_m2)]);
                        
                        % ============================ only continue if have data
                        % for both train and base
                        if N1<Nmin | N2< Nmin
                            %                             disp('SKUIPPED - low N for eitehr base or train');
                            continue
                        end
                        
                        % --- increment motif pair
                        mcount = mcount+1;
                        
                        
                        
                        
                        
                        if getdprime==1
                            % ################################################# DPRIME
                            % ============================= COLLECT BASELINE AND TRAINING FR DIFF
                            % METRICS
                            frmat1 = [motifstats.neurons(nn).motif(mm).SegmentsExtract.FRsmooth_rate_CommonTrialDur];
                            frmat2 = [motifstats.neurons(nn).motif(mmm).SegmentsExtract.FRsmooth_rate_CommonTrialDur];
                            
                            xbins = motifstats.neurons(nn).motif(mm).SegmentsExtract(1).FRsmooth_xbin_CommonTrialDur;
                            
                            % --- prune time bins so matched between two motifs
                            maxbins = min([size(frmat1,1) size(frmat2,1)]);
                            frmat1 = frmat1(1:maxbins, :);
                            frmat2 = frmat2(1:maxbins, :);
                            xbins = xbins(1:maxbins);
                            
                            % ---- 1) Baseline
                            % - set
                            inds1 = baseInds_m1(end-N1+1:end); % use vals closest to onset of WN
                            inds2 = baseInds_m2(end-N2+1:end);
                            
                            % - run
                            mean1 = mean(frmat1(:, inds1),2);
                            mean2 = mean(frmat2(:, inds2),2);
                            
                            var1 = var(frmat1(:, inds1),0, 2);
                            var2 = var(frmat2(:, inds2),0, 2);
                            
                            dprime_base = abs((mean1-mean2)./sqrt((var1+var2)/2));
                            
                            
                            % ------ 2) Train end
                            % - set
                            inds1 = trainInds_m1(end-N1+1:end); % use vals closest to onset of WN
                            inds2 = trainInds_m2(end-N2+1:end);
                            
                            % - run
                            mean1 = mean(frmat1(:, inds1),2);
                            mean2 = mean(frmat2(:, inds2),2);
                            
                            var1 = var(frmat1(:, inds1),0, 2);
                            var2 = var(frmat2(:, inds2),0, 2);
                            
                            dprime_train = abs((mean1-mean2)./sqrt((var1+var2)/2));
                            
                            
                            % =================================== OUTPUT DPRIME
                            Datstruct.exptcount(ecount).neuronnum(nn).motifpair(mcount).DAT_dprime_base = dprime_base;
                            Datstruct.exptcount(ecount).neuronnum(nn).motifpair(mcount).DAT_dprime_train = dprime_train;
                            Datstruct.exptcount(ecount).neuronnum(nn).motifpair(mcount).DAT_xbins = xbins;
                        end
                        
                        
                        
                        
                        
                        % ################################## RUNNING ANOVA
                        clustnum = motifstats.neurons(nn).clustnum;
                        % --- 1st motif
                        numtrials = length(motifstats.neurons(nn).motif(mm).SegmentsExtract);
                        Yall1 = cell(1,numtrials);
                        for tt = 1:numtrials
                            indstmp = motifstats.neurons(nn).motif(mm).SegmentsExtract(tt).spk_Clust==clustnum;
                            sptimes = motifstats.neurons(nn).motif(mm).SegmentsExtract(tt).spk_Times(indstmp);
                            sptimes = sptimes(sptimes < maxduranova);
                            Yall1{tt} = sptimes;
                        end
                        factorlevels1 = 1*ones(1,numtrials);
                        % ---- 2nd motif
                        numtrials = length(motifstats.neurons(nn).motif(mmm).SegmentsExtract);
                        Yall2 = cell(1,numtrials);
                        for tt = 1:numtrials
                            indstmp = motifstats.neurons(nn).motif(mmm).SegmentsExtract(tt).spk_Clust==clustnum;
                            sptimes = motifstats.neurons(nn).motif(mmm).SegmentsExtract(tt).spk_Times(indstmp);
                            sptimes = sptimes(sptimes < maxduranova);
                            Yall2{tt} = sptimes;
                        end
                        factorlevels2 = 2*ones(1,numtrials);
                        
                        if getdprime==1
                            assert(size(frmat1,2) == length(Yall1), 'dfasdf');
                            assert(size(frmat2, 2) == length(Yall2), 'dafasd');
                        end
                        
                        % ====================== BASELINE - running anova
                        inds1 = baseInds_m1(end-N1+1:end); % use vals closest to onset of WN
                        inds2 = baseInds_m2(end-N2+1:end);
                        %--- combine
                        Yall = [Yall1(inds1) Yall2(inds2)];
                        FactorLevels = [factorlevels1(inds1) factorlevels2(inds2)];
                        
                        [Xall1, OmegaSquared1, ~, ~] = ...
                            lt_neural_RunningAnova(Yall, FactorLevels, windowanova, windshiftanova);
                        if (0)
                            figure; hold on;
                            plot(Xall, OmegaSquared, '-k', Xall, Pval, '-r')
                        end
                        
                        
                        % ====================== TRAINING - running anova
                        inds1 = trainInds_m1(end-N1+1:end); % use vals closest to onset of WN
                        inds2 = trainInds_m2(end-N2+1:end);
                        %--- combine
                        Yall = [Yall1(inds1) Yall2(inds2)];
                        FactorLevels = [factorlevels1(inds1) factorlevels2(inds2)];
                        
                        [Xall2, OmegaSquared2, ~, ~] = ...
                            lt_neural_RunningAnova(Yall, FactorLevels, windowanova, windshiftanova);
                        
                        nmax = min([length(Xall1) length(Xall2)]);
                        
                        
                        % ============ OUTPUT
                        Datstruct.exptcount(ecount).neuronnum(nn).motifpair(mcount).DAT_omega_base = OmegaSquared1(1:nmax);
                        Datstruct.exptcount(ecount).neuronnum(nn).motifpair(mcount).DAT_omega_train = OmegaSquared2(1:nmax);
                        Datstruct.exptcount(ecount).neuronnum(nn).motifpair(mcount).DAT_omega_change = OmegaSquared2(1:nmax) - OmegaSquared1(1:nmax);
                        Datstruct.exptcount(ecount).neuronnum(nn).motifpair(mcount).DAT_anovaxbins = Xall2(1:nmax);
                        
                        
                        
                        
                        % ############################################## GET STATS
                        % ABOUT THIS PAIR OF MOTIFS
                        singlesyl_1 = syl1(findstr(syl1, '(')+1);
                        istargcont_1 = SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).STATS_motifsyl(mm).Isthistargcontext;
                        istargsyl_1 = SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).STATS_motifsyl(mm).Isthistargsyl;
                        learndir_1 =  SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).STATS_motifsyl(mm).Thislearndir;
                        
                        
                        singlesyl_2 = syl2(findstr(syl2, '(')+1);
                        istargcont_2 = SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).STATS_motifsyl(mmm).Isthistargcontext;
                        istargsyl_2 = SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).STATS_motifsyl(mmm).Isthistargsyl;
                        learndir_2 =  SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).STATS_motifsyl(mmm).Thislearndir;
                        
                        
                        if strcmp(singlesyl_1, singlesyl_2) ==1
                            samesyl =1;
                        else
                            samesyl = 0;
                        end
                        istargcontext = [istargcont_1 istargcont_2];
                        istargsyl = [istargsyl_1 istargsyl_2];
                        learndir = [learndir_1 learndir_2];
                        
                        Datstruct.exptcount(ecount).neuronnum(nn).motifpair(mcount).motifs = {syl1 syl2};
                        Datstruct.exptcount(ecount).neuronnum(nn).motifpair(mcount).singlesyls = {singlesyl_1 singlesyl_2};
                        Datstruct.exptcount(ecount).neuronnum(nn).motifpair(mcount).issamesyl = ...
                            samesyl;
                        Datstruct.exptcount(ecount).neuronnum(nn).motifpair(mcount).istargcontext = ...
                            istargcontext;
                        Datstruct.exptcount(ecount).neuronnum(nn).motifpair(mcount).istargsyl = ...
                            istargsyl;
                        Datstruct.exptcount(ecount).neuronnum(nn).motifpair(mcount).learndir = ...
                            learndir;
                        
                    end
                end
                
                % =================== TROUBLESHOOT - compares numebr of
                % extractee pairs with sampel sizes - use to dertermine
                % what to get more labels of.
                if (0)
                    if mcount<50
                        
                        
                        functmp = @(X) sum(X);
                        baseN = cellfun(functmp ,{SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).neuron(nn).DATA.motif.baseInds});
                        trainN = cellfun(functmp, {SwitchStruct.bird(bb).exptnum(ee).switchlist(ss).neuron(nn).DATA.motif.trainInds})
                        
                        lt_figure; hold on;
                        subplot(211); hold on;
                        lt_plot_histogram(baseN);
                        title([birdname '-' exptname '-sw' num2str(ss)]);
                        xlabel('base');
                        subplot(212); hold on;
                        lt_plot_histogram(trainN);
                        xlabel('train');
                        title(['numpairs: ' num2str(mcount)])
                    end
                end
            end
        end
    end
end



