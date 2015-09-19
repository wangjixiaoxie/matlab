%% LT 2/26/15 - Perform bootstrap resampling to calculate confidence intervals of statistics
function BootStats=lt_bootstrap(X,Qfunc,Nsamps);

% X: 1D array of data points
% Qfunc: what function of data to get confidence intervals of? Examples:
%     Qfunc='std'; % to get confidence intervals of standard devation
%     Qfunc='mean'; % conf intervals of mean
%     Qfunc='cv'; % conf intervals of cv
% Nsamps: how many resamples to perform? good value is around 1000;
%     leave blank to get default: 1000;\


%% PARAMS
alpha=0.05; % confidence intervals [alpha/2, 1-alpha/2];
N=length(X);

if ~exist('Nsamps','var');
    Nsamps=1000;
    disp('Using default Nsamps of 1000');
end



%% CALCULATE

    disp(['resampling for ' Qfunc '...'])
% 
% % CI of standard deviation
% if strcmp(Qfunc,'std')==1;
%     Qboot=[]; % will fill with resampled STDs
%     for i=1:Nsamps;
%         Xboot=X(randi(N, N, 1)); % take random samples from X, with replacement, and get same sample size
%         Qboot(i)=std(Xboot);
%     end
%     
%     BootStats.STD=std(Qboot);
%     BootStats.MEAN=mean(Qboot);
%     BootStats.MEDIAN=median(Qboot);
%     BootStats.CI=prctile(Qboot,100*[alpha/2, 1-alpha/2]);
%     
% end

% CI of standard deviation
Qboot=[]; % will fill with resampled STDs
for i=1:Nsamps;
    Xboot=X(randi(N, N, 1)); % take random samples from X, with replacement, and get same sample size
    
    
    if strcmp(Qfunc,'std')==1;
        Qboot(i)=std(Xboot);
    elseif strcmp(Qfunc,'cv')==1;
        Qboot(i)=std(Xboot)/mean(Xboot); % get a CV value for this resampled dataset
    elseif strcmp(Qfunc,'mean')==1
        Qboot(i)=mean(Xboot);
    end
    
end
BootStats.STD=std(Qboot);
BootStats.MEAN=mean(Qboot);
BootStats.MEDIAN=median(Qboot);
BootStats.CI=prctile(Qboot,100*[alpha/2, 1-alpha/2]);


%
% if strcmp(Qfunc,'cv')==1;
%         Qboot=[]; % will fill with resampled values
%     for i=1:Nsamps;
%         Xboot=X(randi(N, N, 1)); % take random samples from X, with replacement, and get same sample size
%         Qboot(i)=std(Xboot)/mean(Xboot); % get a CV value for this resampled dataset
%     end
%     
%     BootStats.STD=std(Qboot);
%     BootStats.MEAN=mean(Qboot);
%     BootStats.MEDIAN=median(Qboot);
%     BootStats.CI=prctile(Qboot,100*[alpha/2, 1-alpha/2]);
% end
%     
%     
% 


%% end

% disp('BootStats:');
% disp(BootStats);



