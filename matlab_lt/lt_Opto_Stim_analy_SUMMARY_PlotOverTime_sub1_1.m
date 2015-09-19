%% SET UP shuffling.
tic
Ncycles=5000; % how many cycles.

%% RUN stuff that will not change with shuffling
for i=1:length(TimeFieldsOfInterest);
    NumEpochs=length(StimEpochs_aligned.timewindow{i}.epoch);
    
    for j=1:NumEpochs;        
        if isfield(StimEpochs_aligned.timewindow{i}.epoch{j}.data,'All_preceding');
            
            % COLLECT STATS on this epoch
            % 1) preceding (nonstim)
            nrends=length(StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.([statfield '_fudgeDST'])); % num renditions
            if ~isempty(NumTrials_otherstats);
                % if there are not enough data, give use warning
                if nrends<NumTrials_otherstats;
                    disp(['warning, not enough trials for preceding (nonstim) for day: ' num2str(ii), ' epoch: ' num2str(k) ' have ' num2str(length(Inds)) ' trials. Will take all available rends']);
                    
                    Yvals=StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.([statfield '_fudgeDST']);
                    Tvals=StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.timevals_fudgeDST_minust0;
                    
                else % have enough trials
                    
                    Yvals=StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.([statfield '_fudgeDST'])(end-NumTrials_otherstats+1:end);
                    Tvals=StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.timevals_fudgeDST_minust0((end-NumTrials_otherstats+1:end));
                    
                end
                
            end
            
            StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.STATS_PossiblySubset.([statfield '_fudgeDST']).mean=mean(Yvals);
            
            % slope
            % by time
            X=[ones(length(Tvals),1) Tvals'];
            Y=Yvals';
            b=regress(Y,X);
            
            StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.STATS_PossiblySubset.([statfield '_fudgeDST']).slope=b(2);
            
            % slope by rendition
            X=[ones(1,length(Tvals)); 1:length(Tvals)]';
            Y=Yvals';
            b=regress(Y,X);
            
            StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.STATS_PossiblySubset.([statfield '_fudgeDST']).slope_by_rend=b(2);
        end
    end
end

%% BELOW IS STUFF THAT WILL BE SHUFFLED
for nn=1:Ncycles;
    
    for i=1:length(TimeFieldsOfInterest);
        
        NumEpochs=length(StimEpochs_aligned.timewindow{i}.epoch);
        
        for j=1:NumEpochs;
            
            if isfield(StimEpochs_aligned.timewindow{i}.epoch{j}.data,'All_preceding');
                
                % 2) Stim epoch (catch and not catch)
                tmpfields={'StimCatch','StimNotCatch'};
                Yvals={};
                Tvals={};
                for kkk=1:length(tmpfields);
                    % get raw values for both epochs, then shuffle them.
                    tmpf=tmpfields{kkk};
                    n=length(StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).([statfield '_fudgeDST'])); % num renditions
                    if ~isempty(NumTrials_reversion);
                        if n<NumTrials_reversion;
                            disp(['warning, not enough trials for ' tmpf ' for day: ' num2str(ii), ' epoch: ' num2str(k) ' have ' num2str(length(Inds)) ' trials. Will continue']);
                            
                            Yvals{kkk}=StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).([statfield '_fudgeDST']);
                            Tvals{kkk}=StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).timevals_fudgeDST_minust0;
                            
                        else % have enough trials
                            
                            Yvals{kkk}=StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).([statfield '_fudgeDST'])(1:NumTrials_reversion);
                            Tvals{kkk}=StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).timevals_fudgeDST_minust0(1:NumTrials_reversion);
                        end
                    end
                    
                end
                
                % SHUFFLE THE VALS
                N1=length(Yvals{1});
                N2=length(Yvals{2});
                
                Inds_shuffle=randperm(N1+N2); % get random indices
                
                % shuffle values
                Ytot=[Yvals{1} Yvals{2}];
                Ytot_shuffle=Ytot(Inds_shuffle);
                
                Ttot=[Tvals{1} Tvals{2}];
                Ttot_shuffle=Ttot(Inds_shuffle);
                
                % redistribute into trial types
                Yvals_shuffle{1}=Ytot_shuffle(1:N1);
                Yvals_shuffle{2}=Ytot_shuffle(N1+1:end);
                
                Tvals_shuffle{1}=Ttot_shuffle(1:N1);
                Tvals_shuffle{2}=Ttot_shuffle(N1+1:end);
                
                % sort in temporal order
                for kkk=1:2;
                    [Tvals_shuffle{kkk}, Inds]=sort(Tvals_shuffle{kkk});
                    Yvals_shuffle{kkk}=Yvals_shuffle{kkk}(Inds);
                end
                
                
                
                % CONTINUE WITH STATS EXTRACTION, using shuffled data.
                for kkk=1:length(tmpfields);
                    tmpf=tmpfields{kkk};
                    
                    StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).STATS_PossiblySubset_Shuffled.([statfield '_fudgeDST']).mean=mean(Yvals_shuffle{kkk});
                    
                    % slope
                    % by time
                    X=[ones(length(Tvals_shuffle{kkk}),1) Tvals_shuffle{kkk}'];
                    Y=Yvals_shuffle{kkk}';
                    b=regress(Y,X);
                    
                    %                     figure; hold on;
                    %                     plot(X(:,2),Y,'o');
                    %                     plot(xlim,b(1) + b(2).*xlim);
                    
                    StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).STATS_PossiblySubset_Shuffled.([statfield '_fudgeDST']).slope=b(2);
                    
                    % by rendition
                    X=[ones(1,length(Tvals_shuffle{kkk})); 1:length(Tvals_shuffle{kkk})]';
                    Y=Yvals_shuffle{kkk}';
                    b=regress(Y,X);
                    
                    StimEpochs_aligned.timewindow{i}.epoch{j}.data.(tmpf).STATS_PossiblySubset_Shuffled.([statfield '_fudgeDST']).slope_by_rend=b(2);
                    
                end
            end
        end
    end
    
    %% EXTRACT DATA INTO FORMAT FOR PERFORMING REGRESSIONS.
    % first extract data into format: cell array, each cell containing one time
    % field, and that contains matrix where rows are trials and columns go:
    
    % for comparing mean values
    for i=1:length(TimeFieldsOfInterest);
        for j=1:length(StimEpochs_aligned.timewindow{i}.epoch);
            Ymean{i}(j,1)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.STATS_PossiblySubset.([statfield '_fudgeDST']).mean;
            Ymean{i}(j,2)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.StimCatch.STATS_PossiblySubset_Shuffled.([statfield '_fudgeDST']).mean;
            Ymean{i}(j,3)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.StimNotCatch.STATS_PossiblySubset_Shuffled.([statfield '_fudgeDST']).mean;
            
        end
    end
    
    % for slopes
    for i=1:length(TimeFieldsOfInterest);
        for j=1:length(StimEpochs_aligned.timewindow{i}.epoch);
            Yslope{i}(j,1)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.STATS_PossiblySubset.([statfield '_fudgeDST']).slope;
            Yslope{i}(j,2)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.StimCatch.STATS_PossiblySubset_Shuffled.([statfield '_fudgeDST']).mean;
            Yslope{i}(j,3)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.StimNotCatch.STATS_PossiblySubset_Shuffled.([statfield '_fudgeDST']).mean;
            Yslope{i}(j,4)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.All_preceding.STATS_PossiblySubset.([statfield '_fudgeDST']).slope_by_rend;
            Yslope{i}(j,5)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.StimCatch.STATS_PossiblySubset_Shuffled.([statfield '_fudgeDST']).slope_by_rend;
            Yslope{i}(j,6)=StimEpochs_aligned.timewindow{i}.epoch{j}.data.StimCatch.STATS_PossiblySubset_Shuffled.([statfield '_fudgeDST']).slope;
        end
    end
    

    %% REGRESSIONS AND COLLECT P AND R
    
    % 1) Recent learning
    for i=1:length(TimeFieldsOfInterest);
        % X is stim catch minus preceding
        X=diff(Ymean{i},1,2);
        X=X(:,1);
        
        % Y is stim not catch minus stim catch
        Y=diff(Ymean{i},1,2);
        Y=Y(:,2);
        
        % perform regression
        X=[ones(length(X),1) X];
        [~,~,~,~,stats]=regress(Y,X);
        
        % Save stats
        ShuffleStats.timewindow{i}.regressions.recent_learning.p(nn)=stats(3);
        ShuffleStats.timewindow{i}.regressions.recent_learning.r2(nn)=stats(1);

        
    % 2) Current learning (prestim)
    % X is mean of pre-stim
        X=Ymean{i}(:,1);
        
        % Y is stim not catch minus stim catch
        Y=diff(Ymean{i},1,2);
        Y=Y(:,2);
        
        % perform regression
        X=[ones(length(X),1) X];
        [~,~,~,~,stats]=regress(Y,X);
        
        % Save stats
        ShuffleStats.timewindow{i}.regressions.current_learning_prestim.p(nn)=stats(3);
        ShuffleStats.timewindow{i}.regressions.current_learning_prestim.r2(nn)=stats(1);

        
    % 3) Current learning (stim catch)
    % X is mean of stim catch
        X=Ymean{i}(:,2);
        
        % Y is stim not catch minus stim catch
        Y=diff(Ymean{i},1,2);
        Y=Y(:,2);
        
        % perform regression
        X=[ones(length(X),1) X];
        [~,~,~,~,stats]=regress(Y,X);
        
        % Save stats
        ShuffleStats.timewindow{i}.regressions.current_learning_stimcatch.p(nn)=stats(3);
        ShuffleStats.timewindow{i}.regressions.current_learning_stimcatch.r2(nn)=stats(1);

        
    
    % 4) Slope of prestim (time)
    % X is slope of pre-stim
        X=Yslope{i}(:,1);
        
        % Y is reversion (stim not catch minus stim catch)
        Y=diff(Yslope{i},1,2);
        Y=Y(:,2);
        
        % perform regression
        X=[ones(length(X),1) X];
        [~,~,~,~,stats]=regress(Y,X);
        
        % Save stats
        ShuffleStats.timewindow{i}.regressions.slope_time_prestim.p(nn)=stats(3);
        ShuffleStats.timewindow{i}.regressions.slope_time_prestim.r2(nn)=stats(1);

        
    % 5) Slope during stim catch (time)
    % X is slope of stim catch
        X=Yslope{i}(:,6);
        
        % Y is reversion (stim not catch minus stim catch)
        Y=diff(Yslope{i},1,2);
        Y=Y(:,2);
        
        % perform regression
        X=[ones(length(X),1) X];
        [~,~,~,~,stats]=regress(Y,X);
        
        % Save stats
        ShuffleStats.timewindow{i}.regressions.slope_time_stimcatch.p(nn)=stats(3);
        ShuffleStats.timewindow{i}.regressions.slope_time_stimcatch.r2(nn)=stats(1);
    
    
    % 6) Slope during prestim (rend)
        % X is slope of pre-stim
        X=Yslope{i}(:,4);
        
        % Y is reversion (stim not catch minus stim catch)
        Y=diff(Yslope{i},1,2);
        Y=Y(:,2);
        
        
        % perform regression
        X=[ones(length(X),1) X];
        [~,~,~,~,stats]=regress(Y,X);
        
        % Save stats
        ShuffleStats.timewindow{i}.regressions.slope_rend_prestim.p(nn)=stats(3);
        ShuffleStats.timewindow{i}.regressions.slope_rend_prestim.r2(nn)=stats(1);
    
    
    
    % 7) Slope during stimcatch (rend)
        % X is slope of stim catch
        X=Yslope{i}(:,5);
        
        % Y is reversion (stim not catch minus stim catch)
        Y=diff(Yslope{i},1,2);
        Y=Y(:,2);
        
        
        % perform regression
        X=[ones(length(X),1) X];
        [~,~,~,~,stats]=regress(Y,X);
        
        % Save stats
        ShuffleStats.timewindow{i}.regressions.slope_rend_stimcatch.p(nn)=stats(3);
        ShuffleStats.timewindow{i}.regressions.slope_rend_stimcatch.r2(nn)=stats(1);
    end
end
toc