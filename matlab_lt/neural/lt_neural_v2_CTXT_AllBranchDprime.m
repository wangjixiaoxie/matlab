function ALLBRANCH = lt_neural_v2_CTXT_AllBranchDprime(ALLBRANCH, Nmin, Niter, DprimeNegVersion)
%% gets dprime for all branches/neurons
% plots comparison of sample sizes
% also extracts mean fr (for each branch)

%%
% === PARAMS:
% Niter = 3; % shuffles
% Nmin = 3; % min sample size;
% % Nmin = ALLBRANCH.alignpos(1).ParamsFirstIter.minN; % min sample size;
% 
% % DprimeNegVersion = 'Wohl'; % DONT USE THIS - is biased to be large,  because of lower sample size (this splits data in half)
% DprimeNegVersion = 'shuff'; % correct version, shuffles and resplits, maintaining sampel size.


%% RUNS

numalign = length(ALLBRANCH.alignpos);

for i=1:numalign
    
    numbirds = length(ALLBRANCH.alignpos(i).bird);
    Nall_dat = [];
    Nall_pos = [];
    Nall_neg = [];
    for ii = 1:numbirds
        numbranches = length(ALLBRANCH.alignpos(i).bird(ii).branch);
        
        
        for bb = 1:numbranches
            
            numneurons = length(ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron);
            
            for nn=1:numneurons
                
                if isempty(ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).yvals)
                    continue
                end
                
                
                motifpredur = ALLBRANCH.alignpos(i).ParamsFirstIter.motifpredur;
                motifpostdur = ALLBRANCH.alignpos(i).ParamsFirstIter.motifpostdur;
                
                
                
                % ================================== CALCULATE RUNNING
                % D-PRIME AND FR
                
                % ------------- ACTUAL DAT
                dattmp = ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).FR;
                numclasses = length(dattmp.classnum);
                
                indstokeep = round(1:1000*(motifpredur+motifpostdur-0.005));
                %
                dprimethisbranch = [];
                
                for cc=1:numclasses
                    N1= size(dattmp.classnum(cc).FRsmooth_rate_CommonTrialDur,2);
                    
                    % --------------------------- GET FR
                    if N1<Nmin
                        continue
                    end
                    
                    frmean = mean(dattmp.classnum(cc).FRsmooth_rate_CommonTrialDur(indstokeep, :),2);
                    ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).FR.classnum(cc).frmean = frmean;
%                     frRunningStdOfMean = lt_running_stats(frmean, round(stdbinsize*1000), 1);
%                     frRunningStdOfMean = frRunningStdOfMean.STD';
%                     numnan = length(frmean) - length(frRunningStdOfMean);
%                     frRunningStdOfMean = [nan(floor(numnan/2),1); frRunningStdOfMean ...
%                         ; nan(ceil(numnan/2),1)];
                    
%                     ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).FR.classnum(cc).frRunningStdOfMean ...
%                         = frRunningStdOfMean;
%                     ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).FR.classnum(cc).stdbinsize ...
%                         = single(stdbinsize);
                    
                    
                    
                    % -------------------------- DPRIME
                    Nall_dat = [Nall_dat N1];
                    for ccc = cc+1:numclasses
                        N2 = size(dattmp.classnum(ccc).FRsmooth_rate_CommonTrialDur,2);
                        
                        if N1<Nmin || N2<Nmin
                            continue
                        end
                        
                        % ===== method 2
                        alldat = [dattmp.classnum(cc).FRsmooth_rate_CommonTrialDur(indstokeep,:) ...
                            dattmp.classnum(ccc).FRsmooth_rate_CommonTrialDur(indstokeep,:)];
                        
                        class1mean = nanmean(alldat(:, 1:N1),2);
                        class1var = nanvar(alldat(:, 1:N1)');
                        
                        class2mean = nanmean(alldat(:, N1+1:N1+N2),2);
                        class2var = nanvar(alldat(:, N1+1:N1+N2)');
                        
                        if (0)
                            lt_figure; hold on;
                            shadedErrorBar(1:length(class1mean), class1mean, sqrt(class1var), {'Color', 'r'},1 );
                            shadedErrorBar(1:length(class2mean), class1mean, sqrt(class2var), {'Color', 'b'},1 )
                        end
                        
                        
                        dprime = abs((class1mean-class2mean)./sqrt((class1var'+class2var')/2));
                        % from Wohlgemuth, Sober, 2010 (they took abs
                        % value)
                        
                        % ================ OUTPUT FOR THIS BRANCH
                        dprimethisbranch = [dprimethisbranch dprime];
                    end
                end
                
                ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).DprimeAllPairwise = dprimethisbranch;
                
                
                % ------------------ NEG CONTROL
                % =========== version 1 (shuffle across contexts)
                if strcmp(DprimeNegVersion, 'shuff')
                    dattmp = ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).FR;
                    numclasses = length(dattmp.classnum);
                    
                    %                 DprimeAll = [];
                    dprimethisbranch = [];
                    
                    for cc=1:numclasses
                        N1= size(dattmp.classnum(cc).FRsmooth_rate_CommonTrialDur,2);
                        Nall_neg = [Nall_neg N1];
                        if N1<Nmin
                            continue
                        end
                        for ccc = cc+1:numclasses
                            N2 = size(dattmp.classnum(ccc).FRsmooth_rate_CommonTrialDur,2);
                            
                            
                            if N1<Nmin || N2<Nmin
                                continue
                            end
                            
                            % -- actual dat
                            alldat = [dattmp.classnum(cc).FRsmooth_rate_CommonTrialDur(indstokeep,:) ...
                                dattmp.classnum(ccc).FRsmooth_rate_CommonTrialDur(indstokeep,:)];
                            
                            % ---- shuffle
                            DprimeShuffAll = [];
                            for nshuf = 1:Niter
                                indshuff = randperm(N1+N2);
                                %                              indshuff = 1:(N1+N2);
                                alldatshuff = alldat(:, indshuff);
                                
                                class1mean = nanmean(alldatshuff(:, 1:N1),2);
                                class1var = nanvar(alldatshuff(:, 1:N1)');
                                
                                class2mean = nanmean(alldatshuff(:, N1+1:N1+N2),2);
                                class2var = nanvar(alldatshuff(:, N1+1:N1+N2)');
                                
                                %                             % --- cut slightly shorter, so no mismatch
                                %                             class1mean = class1mean(indstokeep);
                                %                             class1var = class1var(indstokeep);
                                %                             class2mean = class2mean(indstokeep);
                                %                             class2var = class2var(indstokeep);
                                %
                                if (0)
                                    lt_figure; hold on;
                                    shadedErrorBar(1:length(class1mean), class1mean, sqrt(class1var), {'Color', 'r'},1 );
                                    shadedErrorBar(1:length(class2mean), class1mean, sqrt(class2var), {'Color', 'b'},1 )
                                end
                                dprime = abs((class1mean-class2mean)./sqrt((class1var'+class2var')/2));
                                % from Wohlgemuth, Sober, 2010 (they took abs
                                % value)
                                
                                DprimeShuffAll = [DprimeShuffAll dprime];
                            end
                            
                            dprimethisbranch = [dprimethisbranch mean(DprimeShuffAll,2)];
                        end
                    end
                elseif strcmp(DprimeNegVersion, 'Wohl')
                    % then within each syl in stereotyped context, shuffle
                    dattmp = ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).FR;
                    %                     dattmp = ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).FR;
                    numclasses = length(dattmp.classnum);
                    
                    %                 DprimeAll = [];
                    dprimethisbranch = [];
                    
                    for cc=1:numclasses
                        N1= floor(size(dattmp.classnum(cc).FRsmooth_rate_CommonTrialDur,2)./2);
                        N2 = ceil(size(dattmp.classnum(cc).FRsmooth_rate_CommonTrialDur,2)./2);
                        
                        Nall_neg = [Nall_neg N1];
                        if N1<Nmin
                            continue
                        end
                        
                        if N1<Nmin || N2<Nmin
                            continue
                        end
                        
                        % -- actual dat
                        alldat = dattmp.classnum(cc).FRsmooth_rate_CommonTrialDur(indstokeep,:);
                        
                        % ---- shuffle
                        DprimeShuffAll = [];
                        for nshuf = 1:Niter
                            indshuff = randperm(N1+N2);
                            %                              indshuff = 1:(N1+N2);
                            alldatshuff = alldat(:, indshuff);
                            
                            class1mean = nanmean(alldatshuff(:, 1:N1),2);
                            class1var = nanvar(alldatshuff(:, 1:N1)');
                            
                            class2mean = nanmean(alldatshuff(:, N1+1:N1+N2),2);
                            class2var = nanvar(alldatshuff(:, N1+1:N1+N2)');
                            
                            %                             % --- cut slightly shorter, so no mismatch
                            %                             class1mean = class1mean(indstokeep);
                            %                             class1var = class1var(indstokeep);
                            %                             class2mean = class2mean(indstokeep);
                            %                             class2var = class2var(indstokeep);
                            %
                            if (0)
                                lt_figure; hold on;
                                shadedErrorBar(1:length(class1mean), class1mean, sqrt(class1var), {'Color', 'r'},1 );
                                shadedErrorBar(1:length(class2mean), class1mean, sqrt(class2var), {'Color', 'b'},1 )
                            end
                            dprime = abs((class1mean-class2mean)./sqrt((class1var'+class2var')/2));
                            % from Wohlgemuth, Sober, 2010 (they took abs
                            % value)
                            
                            DprimeShuffAll = [DprimeShuffAll dprime];
                        end
                        
                        dprimethisbranch = [dprimethisbranch mean(DprimeShuffAll,2)];
                        
                    end
                    
                end
                ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).DprimeAllPairwise_Neg = dprimethisbranch;
                ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).DprimeAllPairwise_Neg_Niter = Niter;
                
                
                % ------------------- POS CONTROL
                dattmp = ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).FR_POSCONTR;
                numclasses = length(dattmp.classnum);
                dprimethisbranch = [];
                for cc=1:numclasses
                    N1 = size(dattmp.classnum(cc).FRsmooth_rate_CommonTrialDur,2);
                    Nall_pos = [Nall_pos N1];
                    if N1<Nmin
                        continue
                    end
                    
                    % --------------------------- GET FR
                    frmean = mean(dattmp.classnum(cc).FRsmooth_rate_CommonTrialDur(indstokeep, :),2);
                    ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).FR_POSCONTR.classnum(cc).frmean = frmean;
                    
                    for ccc = cc+1:numclasses
                        N2 = size(dattmp.classnum(ccc).FRsmooth_rate_CommonTrialDur,2);
                        
                        if N1<Nmin || N2<Nmin
                            continue
                        end
                        
                        class1mean = nanmean(dattmp.classnum(cc).FRsmooth_rate_CommonTrialDur,2);
                        class1var = nanvar(dattmp.classnum(cc).FRsmooth_rate_CommonTrialDur');
                        class2mean = nanmean(dattmp.classnum(ccc).FRsmooth_rate_CommonTrialDur,2);
                        class2var = nanvar(dattmp.classnum(ccc).FRsmooth_rate_CommonTrialDur');
                        
                        % --- cut slightly shorter, so no mismatch
                        class1mean = class1mean(indstokeep);
                        class1var = class1var(indstokeep);
                        class2mean = class2mean(indstokeep);
                        class2var = class2var(indstokeep);
                        
                        dprime = abs((class1mean-class2mean)./sqrt((class1var'+class2var')/2));
                        % from Wohlgemuth, Sober, 2010 (they took abs
                        % value)
                        
                        dprimethisbranch = [dprimethisbranch dprime];
                    end
                    
                end
                ALLBRANCH.alignpos(i).bird(ii).branch(bb).neuron(nn).DprimeAllPairwise_Pos = dprimethisbranch;
                
            end
            
        end
        
    end
    
    %% ==================== COMPARE SAMPLE SIZES FOR DAT AND POS
    lt_figure; hold on;
    title('comaprison of sample sizes, should be matched for dprime to be interpretable');
    ylabel('N');
    xlabel('dat, pos, neg');
    lt_plot_MultDist({Nall_dat, Nall_pos, Nall_neg}, [1 2 3], 1);
    
end
