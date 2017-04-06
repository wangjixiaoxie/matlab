function lt_neural_v2_ANALY_VocModel_glm(VOCALSTRUCTall)
%% for each time bin, predicts fr based on sequence information (using GLM)



% %% EXTRACT DATA INTO format for modeling
% % lt_figure; hold on;
% % title('anova rel to syl onset');
% Numbirds = length(VOCALSTRUCTall.birds);
% 
% for z=1:Numbirds
%     Numneurons = length(VOCALSTRUCTall.birds(z).neurons);
%     %     z=1
%     for zz=1:Numneurons
%         VOCALSTRUCT = VOCALSTRUCTall.birds(z).neurons(zz).dat;
%         %         zz=12
%         for i=1:numtimebins
%             
%             % ==========
%             tmp = {VOCALSTRUCT.data_vocalrend.Spk_bincounts};
%             SpkCntMat = cell2mat(tmp');
%             SpkCnts = SpkCntMat(:, i);
%             
%             if logtransform==1
%                 SpkCnts = log10(SpkCnts+1);
%             end
%             
%             Syl = [VOCALSTRUCT.data_vocalrend.syl]';
%             Presyl = [VOCALSTRUCT.data_vocalrend.presyl]';
%             Postsyl = [VOCALSTRUCT.data_vocalrend.postsyl]';
%             
%             % ------ 1) LME model
%             if (0)
%                 modelform = 'SpkCnts ~ 1 + Syl + Presyl';
%                 X = table(SpkCnts, Syl, Presyl);
%                 lme = fitlme(X, modelform);
%             end
%             
%             % ----- 2) ANOVA
%             Group = {};
%             for j=1:length(Groupnames)
%                 Group{j} = eval(Groupnames{j});
%             end
%             if useinteraction ==1;
%                 [p, tabletmp] = anovan(SpkCnts, Group, 'display', 'off', ...
%                     'model', 'interaction');
%             else
%                 [p, tabletmp] = anovan(SpkCnts, Group, 'display', 'off');
%             end
%             
%             numgroups = length(Group)+useinteraction;
%             ss_total = tabletmp{3+numgroups,2};
%             ms_error = tabletmp{2+numgroups,5};
%             
%             Eta2_all = [];
%             Omega2_all = [];
%             
%             for j=1:numgroups
%                 ss_effect = tabletmp{1+j,2};
%                 df_effect = tabletmp{1+j,3};
%                 
%                 numerator = ss_effect - df_effect*ms_error;
%                 denominator = ss_total + ms_error;
%                 omega2 = numerator/denominator;
%                 eta2 = ss_effect/ss_total;
%                 
%                 Eta2_all = [Eta2_all eta2];
%                 Omega2_all = [Omega2_all omega2];
%             end
%             
%             ANOVAOUTPUT.timebin(i).Eta2_all = Eta2_all;
%             ANOVAOUTPUT.timebin(i).Omega2_all = Omega2_all;
%         end
%         ANOVAOUTPUT.global.Groupnames = Groupnames;
%         if useinteraction ==1
%             ANOVAOUTPUT.global.Groupnames = [Groupnames 'interaction'];
%         end
%         VOCALSTRUCTall.birds(z).neurons(zz).anova.output = ANOVAOUTPUT;
%         
%         if (0)
%             % =================== PLOT ANOVA RESULTS
%             tmp = {ANOVAOUTPUT.timebin.Omega2_all};
%             OmegaAcrossBins = cell2mat(tmp');
%             
%             numgroups = length(ANOVAOUTPUT.global.Groupnames);
%             plotcols = lt_make_plot_colors(length(Groupnames), 0,0);
%             for i=1:numgroups
%                 groupname = ANOVAOUTPUT.global.Groupnames{i};
%                 
%                 omegavals = OmegaAcrossBins(:,i);
%                 xvals = VOCALSTRUCTall.global.bincenters;
%                 
%                 plot(xvals, omegavals, '-', 'Color' ,plotcols{i});
%             end
%         end
%     end
% end

%% FOR EACH TIME BIN, PREDICT FR BASED ON SEQUENCE/SYL
% CURRENT STATUS:
% enter specific bird, neuron, and time bin. fits poission model to that
% time bin and plots diagnostics. 
% NOTE: uncomment the iteration in order to iterate over birds and neurons.
% TO DO: can add time bin as a predictor

% Numbirds = length(VOCALSTRUCTall.birds);
% 
% for i=1:Numbirds
%     Numneurons = length(VOCALSTRUCTall.birds(i).neurons);
%     
%     for ii=1:Numneurons
    
        i=1;
        ii=1;
        timebin = 35; % replace with iteration
        
        VocalStruct = VOCALSTRUCTall.birds(i).neurons(ii).dat;

        % ====== EXTRACT FEATURES
            tmp = {VocalStruct.data_vocalrend.Spk_bincounts};
            SpkCntMat = cell2mat(tmp');
            SpkCnts = SpkCntMat(:, timebin);
                        
            Syl = [VocalStruct.data_vocalrend.syl]';
            Presyl = [VocalStruct.data_vocalrend.presyl]';
            Postsyl = [VocalStruct.data_vocalrend.postsyl]';
            
            % ======== PERFORM GLM
            % 1) put predictors into table
            X = table(Syl, Presyl, SpkCnts);
            modelspec = 'SpkCnts ~ Presyl';
            mdl = fitglm(X, modelspec, 'Distribution', 'poisson', ...
                'DispersionFlag', 0);
            
            % ==== DIAGNOSTIC PLOTS
            % 1) plot all residuals - overlay with poisson distrib for that
            % mean
            
            % 2) 
            
            % 
            lt_figure; hold on;
            lt_subplot(4,1,1); hold on;
            mdl.plotResiduals('fitted');
            lt_subplot(3,1,2); hold on;
            mdl.plotResiduals('histogram');
             lt_subplot(3,1,3); hold on;
            mdl.plotResiduals('probability');
            
            if (0)
            % plot residuals by hand
%                 assert(mdl.Residuals.Raw == mdl.Fitted.Response - mdl.Variables.SpkCnts, 'asfasd')
                figure; hold on;
                plot(mdl.Fitted.Response, mdl.Residuals.Raw, 'xb');
                
                % --- if poisson, then plot standardized residuals
            end
            
            if (0) % problem, one of the syls will be used as intercept, so this code doesn't work ATM
            % plot all responses per class
            predictorCategory = 'Syl'; % give a name (e.g. 'Syl');
            
            CatLevels = unique(mdl.Variables.(predictorCategory));
            lt_figure; hold on;
            
            for j=1:length(CatLevels)
               catname = CatLevels(j);
                ind = strfind(mdl.Variables.(predictorCategory)', catname);
                
                % find fit mean and avr
                ind = strcmp(mdl.CoefficientNames, [predictorCategory '_' catname]);
                
                mu_fit = mdl.Coefficients.Estimate(ind);
                SE_fit = mdl.Coefficients.SE(ind);
                tmp = mdl.coefCI;
                CI_fit = tmp(ind, :);
                
                lt_plot(j, mu_fit, {'Errors', SE_fit});
                line([j j], CI_fit, 'LineStyle', '--');
                
            end
            
            lt_plot_zeroline;
            end
            
            % =================== PLOT ALL COEFFICIENTS
            lt_figure; hold on;
            lt_subplot(2,1,1); hold on;
            ylabel('coeff (mean, se, CI)');
            
            xcounter = 1;
            xvals = [];
            for j=1:length(mdl.CoefficientNames)
               
                % find fit mean and avr
                mu_fit = mdl.Coefficients.Estimate(j);
                SE_fit = mdl.Coefficients.SE(j);
                tmp = mdl.coefCI;
                CI_fit = tmp(j, :);
                
                if SE_fit==0
                    continue
                end
                
                lt_plot(xcounter, mu_fit, {'Errors', SE_fit});
                line([xcounter xcounter], CI_fit, 'LineStyle', '-');
                
                xvals = [xvals xcounter];
                xcounter = xcounter +1;
            end
            
            lt_plot_zeroline;
            set(gca, 'XTick', 1:xcounter, 'XTickLabel', mdl.CoefficientNames(xvals));
            
%     end       
% end


%% ========= PREDICT SYL/SEQUENCE IDENTITY GIVEN NEURAL





