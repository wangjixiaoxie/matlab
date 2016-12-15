%% LT 11/19/16 - does ANOVA on running FR, for each time bin gets effect size (e.g. omega squared)
% how much does context affect firing rate?

function [Xall, OmegaSquared, EtaSquared, Pval] = lt_neural_RunningAnova(Yspks, FactorLevels, window, windshift)

% NOTE: performs in time bins up to max spike time across trials (i.e. might not 
% maintain same sample size) - 
% maintains same window size

% Yspks = cell array, with each ind containing spike times in sec
% FactorLevels = array, corersponding to Yspks, for levels (e.g. context)
% window=0.015; sec for window.
% windshift=0.002; % how much to slide window

%% === 
numtrials = length(Yspks); 
xmax=max(cellfun(@max, Yspks(~cellfun(@isempty, Yspks))));

xstarts=0:windshift:xmax-window;
xends=window:windshift:xmax;

if isempty(xstarts) && isempty(xends)
    % then is becuase last spike is earlier than end of first window
    xstarts = [0];
    xends = [window];
end

Xall=[]; % trial x bins
OmegaSquared=[];
EtaSquared = [];
Pval = [];

for j=1:length(xstarts)
    
        % for this window, get count
        wmin=xstarts(j);
        wmax=xends(j);
    
        allSpks_tmp = []; % for this timebin
               
    for i=1:numtrials % num trials
        
        spktimes=Yspks{i};

        y=sum(spktimes>wmin & spktimes<=wmax);
        allSpks_tmp = [allSpks_tmp y];
    end
    
    % ======= PERFORM ANOVA
    assert(length(allSpks_tmp) == length(FactorLevels), 'adfasdcasdc');
    [p, tabletmp] = anova1(allSpks_tmp, FactorLevels, 'off');
    
    ss_effect = tabletmp{2,2};
    ss_total = tabletmp{4,2};
    df_effect = tabletmp{2,3};
    ms_error = tabletmp{3,4};
    
    
    % calculate omega squared
    numerator = ss_effect - df_effect*ms_error;
    denominator = ss_total + ms_error;
    omega2 = numerator/denominator;
        
    
    % calculate eta squared
    eta2 = ss_effect / ss_total;
    
    % --- OUTPUT
    Xall=[Xall; mean([wmin wmax])];
    OmegaSquared = [OmegaSquared omega2];
    EtaSquared = [EtaSquared eta2];
    Pval = [Pval p];
    
end



