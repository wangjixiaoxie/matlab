function lt_neural_BatchSmthXCOV(DATSTRUCT, windowmax, binsize, premotorwind, chan1, ...
    chan2, mm, dirfield, pretime)
%% lt 11/14/17 -

%% ========================== cross correlation of smoothed firing in premotor window
% close all;
% windowmax = 0.04;
% binsize = 0.005;
% premotorwind = [-0.03 0.03]; % use for ch14-21 on 11/12, morning.
% % premotorwind = [-0.8 0.02];
% chan1 = 14;
% chan2 = 21;
% mm = 1;
% dirfield = 'base';
% pretime = 0.1; % sec, from onset

%% ###########################e xtract dat
datsm1 = DATSTRUCT.(dirfield).motifnum(mm).DatAll{chan1};
datsm2 = DATSTRUCT.(dirfield).motifnum(mm).DatAll{chan2};
datraw1 = DATSTRUCT.(dirfield).motifnum(mm).DatAllRaw{chan1};
datraw2 = DATSTRUCT.(dirfield).motifnum(mm).DatAllRaw{chan2};
FFall = DATSTRUCT.(dirfield).motifnum(mm).FF;

% =============== cut off to premotor window
t = DATSTRUCT.(dirfield).motifnum(mm).t;
tmp = (pretime+premotorwind);
indstmp = t>=tmp(1) & t<=tmp(2);
datsm1 = datsm1(:,indstmp);
datsm2 = datsm2(:,indstmp);
datraw1 = datraw1(:,indstmp);
datraw2 = datraw2(:,indstmp);



%% ########################### PLOT OVERLAYED ALL TRIALS
figcount=1;
subplotrows=6;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
numtrials = size(datsm1,1);
for i=1:numtrials
    
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(['ch' num2str(chan1) '[r], ch' num2str(chan2) '[k]']);
    plot(datraw1(i,:), 'r');
    plot(datraw2(i,:), 'k');
    plot(datsm1(i,:), 'r', 'LineWidth', 2);
    plot(datsm2(i,:), 'k', 'LineWidth', 2);
    axis tight
end

% ========================== what data to use?
datmat1 = abs(datraw1);
datmat2 = abs(datraw2);




%% ########################### CALCULATE xcov
figcount=1;
subplotrows=6;
subplotcols=3;
fignums_alreadyused=[];
hfigs=[];
hsplots = [];

% --- for each trial, get cross covariance [DAT]
CCall = [];
numtrials = size(datsm1,1);
for i=1:numtrials
    
    [cc, lags] = xcov(datmat1(i,:), datmat2(i,:), ceil(windowmax/binsize), 'coeff');
    
    CCall = [CCall; cc];
    
    %    plot(lags*binsize, cc, '-k');
end

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('DAT');
plot(lags*binsize, CCall, 'Color', [0.7 0.7 0.7]);
lt_plot(lags*binsize, mean(CCall,1), {'Errors', lt_sem(CCall)});
lt_plot_zeroline;

% ================= [SHUFFLE - lag 1 trial]
CCallSHUFFLE = [];
for i=1:numtrials
    
    d1 = datmat1(i,:);
    if i==numtrials
        d2 = datmat2(1,:);
    else
        d2 = datmat2(i+1,:);
    end
    
    [cc, lags] = xcov(d1, d2, ceil(windowmax/binsize), 'coeff');
    
    CCallSHUFFLE = [CCallSHUFFLE; cc];
    
    %    plot(lags*binsize, cc, '-k');
end

[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
hsplots = [hsplots hsplot];
title('SHIFTED')
plot(lags*binsize, CCallSHUFFLE, 'Color', [0.7 0.7 0.7]);
lt_plot(lags*binsize, mean(CCallSHUFFLE,1), {'Errors', lt_sem(CCallSHUFFLE), ...
    'Color', 'r'});
lt_plot_zeroline;

% ----------
linkaxes(hsplots, 'xy');

%% ============== does CC correlate with FF?
[fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
title('cc vs. FF');
y = mean(CCall,2);
% y = max(CCall, [], 2);
x = FFall;
lt_regress(y, x, 1, 0, 1, 1, 'k', 0);