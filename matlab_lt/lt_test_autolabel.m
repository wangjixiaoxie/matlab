clear all; close all;
lt_autolabel

%%
batch='batch.rec_FB';
note='b';
[fnames, sylnum]=lt_jc_chcklbl(batch,note, 0.025,0.025,'','','');

[vlsorfn, vlsorind]=jc_vlsorfn(batch,note,'','');

%% Then check output of above function using evsonganaly

%%
labels_all=[];
for i=1:length(fnames);
    load([fnames{i} '.not.mat']);
    labels_all=[labels_all labels];
end
Inds=find(labels_all~='-');
sylnum_handlabeled=length(labels_all);

if sylnum~=sylnum_handlabeled;
    disp('WARNING - hand labeled # of syl does not match actual number, can cause autolabeling frameshift error - check evsonganaly thresholding');
end

lt_jc_fndlbl(vlsorind,vlsorfn,Inds,labels_all,1)

%% DEBUGGING
% 
% !rm pu11wh87_200814_065926.8340.cbin.not.mat
% !rm pu11wh87_200814_070102.8361.cbin.not.mat
% 
% 
% load pu11wh87_200814_065926.8340.cbin.not.mat
% ONSETS.auto{1}=onsets;
% OFFSETS.auto{1}=offsets;
% 
% load pu11wh87_200814_070102.8361.cbin.not.mat
% ONSETS.auto{2}=onsets;
% OFFSETS.auto{2}=offsets;
% 
% 
% load pu11wh87_200814_065926.8340.cbin.not.mat
% ONSETS.hand{1}=onsets;
% OFFSETS.hand{1}=offsets;
% 
% load pu11wh87_200814_070102.8361.cbin.not.mat
% ONSETS.hand{2}=onsets;
% OFFSETS.hand{2}=offsets;
% 
% %%
% 
% [ONSETS.auto{1}- ONSETS.hand{1}]
% 
% [OFFSETS.auto{1}- OFFSETS.hand{1}]
% 
% [ONSETS.auto{2}(1:14)-ONSETS.hand{2}]
% 
% [OFFSETS.auto{2} OFFSETS.hand{2}]
% 
% figure; hold on;
% plot(ONSETS.auto{2},'r');
% plot(ONSETS.hand{2},'b');
% 
% 
% ONSETS.auto{2}(13:15)
% ONSETS.hand{2}(13:15)
% 
% 

