%%% Calculating FF & 75% percentile of FF from labeled, triggered files (10/04/11)

%% labeling files for calculating FF 
cd /canary/SY/bk90bk46/060612_thrshtst
dirf('*.cbin','batch');
findcatch('batch');
edit batch.catch % 10/04/11 after 10am data only 

% to find out how many files in batch
% this batch contains a file where WN was not given even it was triggered
% by criteria. To see if the triggered syllable was really below the
% threshold (=2775Hz), I'm checking FF for this. 
%findcatch('batch') % since there were 200 files in batch, choose catch f2iles only (0.1 fraction) for labeling
evsonganaly % label the syllable of interest in catch files 

%% Calculating FF among labeled, triggered files
vals_A0606=evtaf_freq('batch.catch_0606', [4500 5500], 'a', 128, 'obs0', 1); 
vals_B0606=evtaf_freq('batch.catch_0606', [4500 5500], 'b', 128, 'obs0', 1); 

% this only detects triggered. so cannot use for catch during WN
% vals=evtaf_freq(bt,fbins,NT,NFFT,CS,USEX);
%
%  bt - batch filefigure; hold on; plot(mean(pc_bk27bk7_baseUD'),'k')
%  fbins - [Min Freq, Max Freq] to search for peaks (in Hz); [2000 3000],
%  lower & longer; [2800 3800], higher & shorter
%  
%  NT - target note, if you labe the syllable of interest as 'a', use 'a'
%  NFFT - length of template 
%  CS - chan spec
%  USEX - if == 1 look in X.rec for trigger times
% returns vals - 
%   vals=[datenum of file , FREQ vals , Note Index , file number]-> not for
%   this file, FREQ vals should be in 6th column.... (SY)
% usage example:
%  vals=evtaf_freq('batch.train',[5000,6000],'a',128,'obs0',0);

% in this function FFs are in the 6th column.1st function is time.
figure; hold on; plot (vals_A(:,2));
figure; plot (vals_B(:,2),'r')
%% Calculating mean/std from FF population
mean (vals_A(:,2))
std (vals_A(:,2))

%% Calculating percentile from FF population
prctile(vals_A(:,2), 30) %30 percentbatchile
prctile(vals_A(:,2), 50) 
prctile(vals_A(:,2), 70) 

%% Calculating FF using findwnote2tw.m
dirf('*.cbin.not.mat','batchnotes_bk90bk46_0606') % making batchnotes 10/04 morning
[fvalsstr_0606b]=findwnote2tw('batchnotes_bk90bk46_0606','b','',0.085,[4500 5500],1024,0,'obs0');
% pitchcontour.m to calculate contour of entire pitchplot(mean(pc_bk27bk7_WN'),'r','LineWidth', 2)

[fvalsstr_forpc_0606b]=findwnote2tw('batchnotes_bk90bk46_0606','b','',-0.016,[4500 5500],8000,0,'obs0');
pc_bk90bk46_0606b=jc_pitchcontourFV(fvalsstr_forpc_0606b,1024,1020,1, 4500,5500,[1 2],'obs0');


%jc_pitchmat1024(fvals,1024,1020,1,2000,2700,[1 2],'obs0');
 % batchnotes=batch; a=syllable of interest; 0.08=80ms inward; [2000 3000]
 % peak freq range; 1024 = No. of points to take; 0=not adding notes,
 % obs0=CBIN file (file type), this can be multiple when recording from
 % neurons as weel as recording song. e.g., obs0=song, obs1=neuron1, obs2=
 % neuron2 etc... 
for i=1:length(fvalsstr);FFn(i)=fvalsstr(i).mxvals(2);end
figure;plot(FFn);hold on;plot(vals,'r')
figure;plot(FFn);hold on;plot(vals(:,6),'r')
corrcoef(FFn,vals(:,6)plot(mean(pc_bk27bk7_WN'),'r','LineWidth', 2)

corrcoef(FFn,vals(:,6))
% if the red and blue lines overlaps, calculate FF may be pretty accurate
% (meaning no wrong labeling etc. ??)

% pitch contour all trials
figure;hold on; plot ((pc_bk90bk46_0606),'k');
plot ((pc_bk90bk46_0608WN),'r');

figure;hold on; plot ((pc_bk90bk46_0606b),'k');
plot ((pc_bk90bk46_0608WNb),'b');

% comparing mean +/- STD of pitch contours before (blue) and during 
figure; hold on; plot(mean(pc_bk90bk46_0606'),'k','LineWidth', 2)
plot(mean(pc_bk90bk46_0606')+std(pc_bk90bk46_0606'), 'k');
plot(mean(pc_bk90bk46_0606')-std(pc_bk90bk46_0606'), 'k');

plot(mean(pc_bk90bk46_0608WN'),'r','LineWidth', 2)
plot(mean(pc_bk90bk46_0608WN')+std(pc_bk90bk46_0608WN'), 'r');
plot(mean(pc_bk90bk46_0608WN')-std(pc_bk90bk46_0608WN'), 'r');

figure; hold on; plot(mean(pc_bk90bk46_0606b'),'k','LineWidth', 2)
plot(mean(pc_bk90bk46_0606b')+std(pc_bk90bk46_0606b'), 'k');
plot(mean(pc_bk90bk46_0606b')-std(pc_bk90bk46_0606b'), 'k');

plot(mean(pc_bk90bk46_0608WNb'),'b','LineWidth', 2)figure;hold on;
plot((timing4(fvalsstr_0606b)),(mean(pc_bk90bk46_0606b(180:220,:))),'k.');
plot (3618, mean(mean(pc_bk90bk46_0606b(180:220,:))), 'ko');
plot (3618, mean(mean(pc_bk90bk46_0606b(180:220,:)))-std(mean(pc_bk90bk46_0606b(180:220,:))), 'k-');
plot (3618, mean(mean(pc_bk90bk46_0606b(180:220,:)))+std(mean(pc_bk90bk46_0606b(180:220,:))), 'k-');

plot((timing4(fvalsstr_0608WNb)),(mean(pc_bk90bk46_0608WNb(180:220,:))),'b.');
plot (3618, mean(mean(pc_bk90bk46_0608WNb(180:220,:))), 'bo');
plot (3618, mean(mean(pc_bk90bk46_0608WNb(180:220,:)))-std(mean(pc_bk90bk46_0608WNb(180:220,:))), 'b-');
plot (3618, mean(mean(pc_bk90bk46_0608WNb(180:220,:)))+std(mean(pc_bk90bk46_0608WNb(180:220,:))), 'b-');
plot(mean(pc_bk90bk46_0608WNb')+std(pc_bk90bk46_0608WNb'), 'b');
plot(mean(pc_bk90bk46_0608WNb')-std(pc_bk90bk46_0608WNb'), 'b');

%% FF plot
% syllable a
figure;hold on;
plot((timing4(fvalsstr_0606)),(mean(pc_bk90bk46_0606(180:220,:))),'k.');
plot (3618, mean(mean(pc_bk90bk46_0606(180:220,:))), 'ko');
plot (3618, mean(mean(pc_bk90bk46_0606(180:220,:)))-std(mean(pc_bk90bk46_0606(180:220,:))), 'k-');
plot (3618, mean(mean(pc_bk90bk46_0606(180:220,:)))+std(mean(pc_bk90bk46_0606(180:220,:))), 'k-');

plot((timing4(fvalsstr_0608WN)),(mean(pc_bk90bk46_0608WN(180:220,:))),'r.');
plot (3618, mean(mean(pc_bk90bk46_0608WN(180:220,:))), 'ro');
plot (3618, mean(mean(pc_bk90bk46_0608WN(180:220,:)))-std(mean(pc_bk90bk46_0608WN(180:220,:))), 'r-');
plot (3618, mean(mean(pc_bk90bk46_0608WN(180:220,:)))+std(mean(pc_bk90bk46_0608WN(180:220,:))), 'r-');

% syllable b
figure;hold on;
plot((timing4(fvalsstr_0606b)),(mean(pc_bk90bk46_0606b(180:220,:))),'k.');
plot (3618, mean(mean(pc_bk90bk46_0606b(180:220,:))), 'ko');
plot (3618, mean(mean(pc_bk90bk46_0606b(180:220,:)))-std(mean(pc_bk90bk46_0606b(180:220,:))), 'k-');
plot (3618, mean(mean(pc_bk90bk46_0606b(180:220,:)))+std(mean(pc_bk90bk46_0606b(180:220,:))), 'k-');

plot((timing4(fvalsstr_0608WNb)),(mean(pc_bk90bk46_0608WNb(180:220,:))),'b.');
plot (3618, mean(mean(pc_bk90bk46_0608WNb(180:220,:))), 'bo');
plot (3618, mean(mean(pc_bk90bk46_0608WNb(180:220,:)))-std(mean(pc_bk90bk46_0608WNb(180:220,:))), 'b-');
plot (3618, mean(mean(pc_bk90bk46_0608WNb(180:220,:)))+std(mean(pc_bk90bk46_0608WNb(180:220,:))), 'b-');


%Variability
figure;hold on;plot(prctile(pc_r4w53b_0530UD',84)-prctile(pc_r4w53b_0530UD',16),'k');
plot(prctile(pc_r4w53b_0530FD',84)-prctile(pc_r4w53b_0530FD',16),'r');
plot(prctile(pc_r4w53_0530am',84)-prctile(pc_r4w53_0530am',16),'b');
