function [Ypredicted, ConfMat, Accuracy, Sensitivity, PiYActual] = lt_neural_v2_QUICK_classify(Xinput, ...
    Yinput, method, rebalance, imbalance_thr, adasyn_beta)
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

%% change, so that categories ordered 1,2, ... (will convert back at end)

[Yinput_G, GN, GL] = grp2idx(Yinput);


%%

% == performs leave one out cross validation in making prediction


Ypredicted = [];
Pihat_all = [];

for xx =1:size(Xinput,1)
    
    X_train = Xinput([1:xx-1 xx+1:end], :);
    Y_train = Yinput_G([1:xx-1 xx+1:end]);
    
    X_test = Xinput(xx, :);
    
    % --- make sure is full rank
    if rank(X_train)< size(X_train,2)
        Ypredicted = [];
        Pihat_all = [];
        disp('NOT FULL RANK - failed');
        break
    end
    
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
               
               [x, y] = ...
                   ADASYN(X_train, Ytmp, betatmp);

               Xsyn = [Xsyn; x];
               Ysyn = [Ysyn; nn*ones(size(x,1),1)];
           end
           
           % ---- add to actual data
           X_train = [X_train; Xsyn];
           Y_train = [Y_train; Ysyn];
       end
    end
    
    % - training
    % --- matlab version (no regularization)
    if strcmp(method, 'matlab')
        
        [B, dev, stats] = mnrfit(X_train, Y_train);
        
        % - test
        pihat = mnrval(B, X_test);
        [~, ypred] = max(pihat);
        
        
    elseif strcmp(method, 'glmnet');
        
        %         cvfit = cvglmnet(X_train, Y_train, 'multinomial', options);
        fit = glmnet(X_train, Y_train, 'multinomial', options);
        
        %         if rand<0.5
        %             cvglmnetPlot(cvfit)
        %             pause
        %         end
        
        %         ypred = glmnetPredict(fit, X_test, options.lambda, 'class');
        pihat = glmnetPredict(fit, X_test, options.lambda, 'response');
        [~, ypred] = max(pihat);
    end
    
    Pihat_all = [Pihat_all; pihat];
    Ypredicted = [Ypredicted; ypred];
end



    % ================ convert back to original categorical matrices
    Ypredicted = GL(Ypredicted);


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
