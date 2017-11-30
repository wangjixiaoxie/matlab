function [Ypredicted, ConfMat, Accuracy, Sensitivity, PiYActual] = lt_neural_v2_QUICK_classify(Xinput, ...
    Yinput, method, rebalance, imbalance_thr, adasyn_beta, numfold)
%% lt 9/26/17 - changed from LOO to 10-fold (was getting too slow)

% numfold = 10;

%% lt 8/11/17 - multinomial classifier using firing rate vector as input
% X, firing rate matrix (trials x tbins)
% Y, categorical variable (using integers, but not ordered) (trials x 1)

if ~exist('method', 'var')
    method = 'glmnet'; % includes regulaization
end

if strcmp(method, 'glmnet')
    options = glmnetSet;
    options.alpha = 0.75; % 1 is L1, 0 is L2 (in between is mixing)
    options.lambda = exp(-3); % based on empirically looking at various window sizes (40ms to 100ms) and binsizes (2ms to 20ms)
end

if ~exist('rebalance', 'var')
    rebalance = 0; % % if 1, then will upsample minority classes on each run of classification (if imbalanced)
    
end

if ~exist('imbalance_thr', 'var')
    imbalance_thr=0.5; % if 0.5, then rebalances if smallest class if <50% sample size of largest class.
end


if ~exist('adasyn_beta', 'var')
    adasyn_beta=0.8;
end
%     adasyn_beta = 0.6; % fraction to reduce gap between minor and major calsses (if doing rebalancing)


Xinput = double(Xinput);
%% change, so that categories ordered 1,2, ... (will convert back at end)

[Yinput_G, GN, GL] = grp2idx(Yinput);


%% =================== get inds for cross validation

if (0)
    % method 1) do separtely for each class, so no change that one iteration is
    % missing a class entirely.
    tmp = tabulate(Yinput);
    tmp = cell2mat(tmp(:,2));
    
    kfoldinds = [];
    for i=1:length(tmp)
        indstmp = crossvalind('Kfold', tmp(i), numfold);
        
        kfoldinds = [kfoldinds; indstmp];
    end
    
    assert(all(size(kfoldinds) == [length(Yinput) 1]), 'asfsda');
    
    % method 2) do on all data at once
    kfoldinds = crossvalind('Kfold', length(Yinput), numfold);
end

% -- current method, makes sure at least one of each group is in training
% set.

% tmp = cvpartition(Yinput, 'kfold', numfold)


%%
if (0) % OLD VERSION - doesn't work on matlab 2013, so changed to other version using cvpartition
    % NOTE!!!!!!!!!!!: if use old version then need to modify by puitting
    % if,else, based on var of predictors.
    % == performs leave one out cross validation in making prediction
    
kfoldinds = crossvalind('Kfold', Yinput, numfold);
    
    Ypredicted = nan(length(Yinput),1);
    Pihat_all = nan(length(Yinput),length(unique(Yinput)));
    
    tmpsum = 0;
    
    kfoldlist = unique(kfoldinds);
    
    for xx = kfoldlist'
        X_test = Xinput(kfoldinds==xx, :);
        
        X_train = Xinput(kfoldinds~=xx, :);
        Y_train = Yinput_G(kfoldinds~=xx);
        
        %     xx =1:size(Xinput,1)
        %
        %     X_train = Xinput([1:xx-1 xx+1:end], :);
        %     Y_train = Yinput_G([1:xx-1 xx+1:end]);
        %
        %     X_test = Xinput(xx, :);
        %
        
        
        
        % --- make sure is full rank
        %         if rank(X_train)< size(X_train,2)
        %             Ypredicted = [];
        %             Pihat_all = [];
        %             disp('NOT FULL RANK - failed');
        %             break
        %         end
        
        
        % --- if imbalanced, make balanced
        if rebalance==1
            % -- determine how unbalanced this is
            Nall = grpstats(Y_train, Y_train, 'numel');
            if min(Nall)/max(Nall) < imbalance_thr
                
                % --- rebalance by generating syntehtic data using ADASYN (puts
                % more data at boundaries
                
                % for each minority class, generate synthetic data
                Xsyn = [];
                Ysyn = [];
                for nn=1:length(Nall)
                    
                    if Nall(nn) == max(Nall)
                        % skip, this is max class
                        continue
                    end
                    
                    indsMinor = Y_train == nn;
                    indsMajor = ~indsMinor;
                    
                    % recode, so that minor class is 0 and major class is 1
                    Ytmp = Y_train;
                    Ytmp(indsMinor) = 0;
                    Ytmp(indsMajor) = 1;
                    
                    % figure out how many new samples needed
                    numsampsneeded = max(Nall) - Nall(nn);
                    
                    betatmp = adasyn_beta * numsampsneeded/(sum(indsMajor)-sum(indsMinor)); % downscale beta, since majority class is all other classes
                    
                    % if all
                    if all(std(X_train(indsMinor,:),1) ==0)
                        % then no variability ... do dumb think just replicate
                        % all points
                        x = repmat(X_train(find(indsMinor,1,'first'),:), floor(adasyn_beta * numsampsneeded),1);
                        disp('SYNTHESIZING NEW POINTS BY COPYING (SINCE NO VARIABILITY)');
                    else
                        % syntheisze new datapoints
                        [x, ~] = ...
                            ADASYN(X_train, Ytmp, betatmp);
                    end
                    
                    Xsyn = [Xsyn; x];
                    Ysyn = [Ysyn; nn*ones(size(x,1),1)];
                end
                
                % ---- add to actual data
                X_train = [X_train; Xsyn];
                Y_train = [Y_train; Ysyn];
            end
        end
        
        
        
        % ########################################### training
        % --- matlab version (no regularization)
        if strcmp(method, 'matlab')
            
            [B, dev, stats] = mnrfit(double(X_train), Y_train);
            
            % - test
            pihat = mnrval(B, X_test);
            [~, ypred] = max(pihat);
            
            
        elseif strcmp(method, 'glmnet')
            
            %         cvfit = cvglmnet(X_train, Y_train, 'multinomial', options);
            
            % === version that does try
            %             donedone=0; % to catch error and retry until works.
            %             while donedone==0
            %                 try
            %                     fit = glmnet(X_train, Y_train, 'multinomial', options);
            %                     donedone=1;
            %                 catch err
            %                     donedone=0;
            %                     disp('-- RETRYING glmnet, ran into error...');
            %                 end
            %             end
            
            % ===
            fit = glmnet(X_train, Y_train, 'multinomial', options);
            
            %         if rand<0.5
            %             cvglmnetPlot(cvfit)
            %             pause
            %         end
            
            %         ypred = glmnetPredict(fit, X_test, options.lambda, 'class');
            pihat = glmnetPredict(fit, X_test, options.lambda, 'response');
            [~, ypred] = max(pihat, [], 2);
        end
        
        %     Pihat_all = [Pihat_all; pihat];
        %     Ypredicted = [Ypredicted; ypred];
        Pihat_all(kfoldinds==xx,:) = pihat;
        Ypredicted(kfoldinds==xx) = ypred;
        
    end
    
    
    
    % ================ convert back to original categorical matrices
    Ypredicted = GL(Ypredicted);
    
    
    
    %%
else
    c = cvpartition(Yinput, 'KFold', numfold);
    
    % ====================
    Ypredicted = nan(length(Yinput),1);
    Pihat_all = nan(length(Yinput),length(unique(Yinput)));
    
    
    for i = 1:c.NumTestSets
        
        X_test = Xinput(test(c,i), :);
        X_train = Xinput(training(c,i), :);
        
        Y_train = Yinput_G(training(c,i));
        
        %     xx =1:size(Xinput,1)
        %
        %     X_train = Xinput([1:xx-1 xx+1:end], :);
        %     Y_train = Yinput_G([1:xx-1 xx+1:end]);
        %
        %     X_test = Xinput(xx, :);
        %
        
        
        
        % --- make sure is full rank
        %         if rank(X_train)< size(X_train,2)
        %             Ypredicted = [];
        %             Pihat_all = [];
        %             disp('NOT FULL RANK - failed');
        %             break
        %         end
        
        
        % --- if imbalanced, make balanced
        if rebalance==1
            % -- determine how unbalanced this is
            Nall = grpstats(Y_train, Y_train, 'numel');
            if min(Nall)/max(Nall) < imbalance_thr
                
                % --- rebalance by generating syntehtic data using ADASYN (puts
                % more data at boundaries
                
                % for each minority class, generate synthetic data
                Xsyn = [];
                Ysyn = [];
                for nn=1:length(Nall)
                    
                    if Nall(nn) == max(Nall)
                        % skip, this is max class
                        continue
                    end
                    
                    indsMinor = Y_train == nn;
                    indsMajor = ~indsMinor;
                    
                    % recode, so that minor class is 0 and major class is 1
                    Ytmp = Y_train;
                    Ytmp(indsMinor) = 0;
                    Ytmp(indsMajor) = 1;
                    
                    % figure out how many new samples needed
                    numsampsneeded = max(Nall) - Nall(nn);
                    
                    betatmp = adasyn_beta * numsampsneeded/(sum(indsMajor)-sum(indsMinor)); % downscale beta, since majority class is all other classes
                    
                    % if all
                    if all(std(X_train(indsMinor,:),1) ==0)
                        % then no variability ... do dumb think just replicate
                        % all points
                        x = repmat(X_train(find(indsMinor,1,'first'),:), floor(adasyn_beta * numsampsneeded),1);
                        disp('SYNTHESIZING NEW POINTS BY COPYING (SINCE NO VARIABILITY)');
                    else
                        % syntheisze new datapoints
                        [x, ~] = ...
                            ADASYN(X_train, Ytmp, betatmp);
                    end
                    
                    Xsyn = [Xsyn; x];
                    Ysyn = [Ysyn; nn*ones(size(x,1),1)];
                end
                
                % ---- add to actual data
                X_train = [X_train; Xsyn];
                Y_train = [Y_train; Ysyn];
            end
        end
        
        
        %% ======= troubleshoot (make sure no group has 0 var for all predictors)
        if (0)
            [grpvar] = grpstats(X_train, Y_train, 'var');
            if all(mean(grpvar==0, 2)==1)
                % then at all groups have 0 var...
                keyboard
            end
            
            
            % --- any columns with 0 var?
            if all(var(X_train)==0)
                keyboard
            end
        end
        
        %%
        
        
        % ########################################### training
        % --- matlab version (no regularization)
        if strcmp(method, 'matlab')
            
            [B, dev, stats] = mnrfit(double(X_train), Y_train);
            
            % - test
            pihat = mnrval(B, X_test);
            [~, ypred] = max(pihat);
            
            
        elseif strcmp(method, 'glmnet')
            
            if all(var(X_train)==0)
                % then bad, since no variance for any predictor. glmnet
                % will crash. solution: skip glmnet and predict using group
                % frequencies
                
                tmp = tabulate(Y_train);
                pihat = repmat(tmp(:,3)'./100, size(X_test,1), 1); % frequencies of each group
                
            else
                % === version that does try
                donedone=0; % to catch error and retry until works.
                while donedone==0
                    try
                        fit = glmnet(X_train, Y_train, 'multinomial', options);
                        donedone=1;
                    catch err
                        donedone=0;
                        disp('-- RETRYING glmnet, ran into error...');
                        keyboard
                    end
                end
                
                % ===
                %                     fit = glmnet(X_train, Y_train, 'multinomial', options);
                
                %         if rand<0.5
                %             cvglmnetPlot(cvfit)
                %             pause
                %         end
                
                %         ypred = glmnetPredict(fit, X_test, options.lambda, 'class');
                pihat = glmnetPredict(fit, X_test, options.lambda, 'response');
            end
            
            [~, ypred] = max(pihat, [], 2);
        end
        
        %     Pihat_all = [Pihat_all; pihat];
        %     Ypredicted = [Ypredicted; ypred];
        Pihat_all(test(c,i),:) = pihat;
        Ypredicted(test(c,i)) = ypred;
        
    end
    
    
    
    % ================ convert back to original categorical matrices
    Ypredicted = GL(Ypredicted);
end

%%
%
% % == performs leave one out cross validation in making prediction
%
%
% Ypredicted = [];
% Pihat_all = [];
%
%
%
% for xx =1:size(Xinput,1)
%
%     X_train = Xinput([1:xx-1 xx+1:end], :);
%     Y_train = Yinput_G([1:xx-1 xx+1:end]);
%
%     X_test = Xinput(xx, :);
%
%     % --- make sure is full rank
%     if rank(X_train)< size(X_train,2)
%         Ypredicted = [];
%         Pihat_all = [];
%         disp('NOT FULL RANK - failed');
%         break
%     end
%
%     % --- if imbalanced, make balanced
%     if rebalance==1
%        % -- determine how unbalanced this is
%        Nall = grpstats(Y_train, Y_train, 'numel');
%        if min(Nall)/max(Nall) < imbalance_thr
%
%            % --- rebalance by generating syntehtic data using ADASYN (puts
%            % more data at boundaries
%
%            % for each minority class, generate synthetic data
%            Xsyn = [];
%            Ysyn = [];
%            for nn=1:length(Nall)
%
%                if Nall(nn) == max(Nall)
%                    % skip, this is max class
%                    continue
%                end
%
%                indsMinor = Y_train == nn;
%                indsMajor = ~indsMinor;
%
%                % recode, so that minor class is 0 and major class is 1
%                Ytmp = Y_train;
%                Ytmp(indsMinor) = 0;
%                Ytmp(indsMajor) = 1;
%
%                % figure out how many new samples needed
%                numsampsneeded = max(Nall) - Nall(nn);
%
%                betatmp = adasyn_beta * numsampsneeded/(sum(indsMajor)-sum(indsMinor)); % downscale beta, since majority class is all other classes
%
%                % if all
%                if all(std(X_train(indsMinor,:),1) ==0)
%                    % then no variability ... do dumb think just replicate
%                    % all points
%                     x = repmat(X_train(find(indsMinor,1,'first'),:), floor(adasyn_beta * numsampsneeded),1);
%                     disp('SYNTHESIZING NEW POINTS BY COPYING (SINCE NO VARIABILITY)');
%                else
%                   % syntheisze new datapoints
%                [x, ~] = ...
%                    ADASYN(X_train, Ytmp, betatmp);
%                end
%
%                Xsyn = [Xsyn; x];
%                Ysyn = [Ysyn; nn*ones(size(x,1),1)];
%            end
%
%            % ---- add to actual data
%            X_train = [X_train; Xsyn];
%            Y_train = [Y_train; Ysyn];
%        end
%     end
%
%     % - training
%     % --- matlab version (no regularization)
%     if strcmp(method, 'matlab')
%
%         [B, dev, stats] = mnrfit(double(X_train), Y_train);
%
%         % - test
%         pihat = mnrval(B, X_test);
%         [~, ypred] = max(pihat);
%
%
%     elseif strcmp(method, 'glmnet')
%
%         %         cvfit = cvglmnet(X_train, Y_train, 'multinomial', options);
%         donedone=0; % to catch error and retry until works.
%         while donedone==0
%             try
%                 fit = glmnet(X_train, Y_train, 'multinomial', options);
%                 donedone=1;
%             catch err
%                 donedone=0;
%                 disp('-- RETRYING glmnet, ran into error...');
%             end
%         end
%         %         if rand<0.5
%         %             cvglmnetPlot(cvfit)
%         %             pause
%         %         end
%
%         %         ypred = glmnetPredict(fit, X_test, options.lambda, 'class');
%         pihat = glmnetPredict(fit, X_test, options.lambda, 'response');
%         [~, ypred] = max(pihat);
%     end
%
%     Pihat_all = [Pihat_all; pihat];
%     Ypredicted = [Ypredicted; ypred];
% end
%
%
%
%     % ================ convert back to original categorical matrices
%     Ypredicted = GL(Ypredicted);


%%
ConfMat = [];
Accuracy = [];
Sensitivity = [];
PiYActual = [];
if ~isempty(Ypredicted)
    
    % ===================== extract prob assigned to the actual response
    PiYActual = Pihat_all(sub2ind(size(Pihat_all), 1:size(Pihat_all,1), Yinput_G')); % prob assigned to the actual outcome
    
    
    % ===================== get stats
    %                    - get confusion matrix
    ConfMat = confusionmat(Yinput, Ypredicted); % actual x predicted
    
    % - classification accuracy metrics
    % accuracy (normalize over all observations)
    N = sum(ConfMat(:));
    Accuracy = trace(ConfMat./N);
    
    % - sensitivity (for actual classes, what probability
    % get correct?)
    tmp = repmat(sum(ConfMat, 2), 1, size(ConfMat, 2));
    Sensitivity = mean(diag(ConfMat./tmp));
end
