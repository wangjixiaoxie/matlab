% below is for wh73pk61

clear all; close all;
% PARAMS
ampThresh=21000;
min_dur=13;
min_int=1;
batch='batch.rec_FB';

templatefile = 'autolabel_templ_b1_v2.dat';
syl.pre='';
syl.post='';
syl.targ='b';

overwrite_notmat=1;


% TEMPLATE SETTINGS
cntrng(1).MIN=1;
cntrng(1).MAX=3;
cntrng(1).NOT=0;
cntrng(1).MODE=1;
cntrng(1).TH=1;
cntrng(1).AND=0;
cntrng(1).BTMIN=0;

cntrng(2).MIN=1;
cntrng(2).MAX=3;
cntrng(2).NOT=0;
cntrng(2).MODE=1;
cntrng(2).TH=2.2;
cntrng(2).AND=0;
cntrng(2).BTMIN=0;

cntrng(3).MIN=1;
cntrng(3).MAX=3;
cntrng(3).NOT=0;
cntrng(3).MODE=1;
cntrng(3).TH=2.2;
cntrng(3).AND=1;
cntrng(3).BTMIN=0;

cntrng(4).MIN=0;
cntrng(4).MAX=10;
cntrng(4).NOT=0;
cntrng(4).MODE=0;
cntrng(4).TH=3;
cntrng(4).AND=0;
cntrng(4).BTMIN=0;

refract=0.2;

% 1) RUN
lt_jc_autoLabel(batch,templatefile,cntrng,syl,refract,ampThresh,min_dur,min_int,overwrite_notmat)


% 2) TO CHECK ACCURACY - isolates all syls labeled and saves to .wav - open and check by eye.
[fnames, sylnum]=lt_jc_chcklbl(batch,syl.targ, 0.025,0.025,'','','');
[vlsorfn vlsorind]=jc_vlsorfn(batch,syl.targ,'','');

% Troubleshooting, enter syl name you want to find song of.
% [fnames, sylnum]=lt_jc_chcklbl(batch,'x', 0.025,0.025,'','','');
% [vlsorfn vlsorind]=jc_vlsorfn(batch,'x','','');
   
disp('DONE!');

    %% Check output using evsonganaly
    
evsonganaly;
    
    
%% After using evsonganaly, run below to replace with corrected things

%     X1=load('syllwv1.wav.not.mat');
%     X2=load('syllwv2.wav.not.mat');
%     XcatLabels=strcat(X1.labels, X2.labels);
%     [Inds, NewSyl]=regexp(XcatLabels,'\w','start','match');
% %     newsyl=strjoin(NewSyl,'\0');
% 
%     lt_jc_fndlbl(vlsorind, vlsorfn, Inds, NewSyl,1)

% collect all labels
labels_all=[];
for i=1:length(fnames);
    load([fnames{i} '.not.mat']);
    labels_all=[labels_all labels];
end

% Any thing not '-' will be replaced by the entered label.
Inds=find(labels_all~='-');
sylnum_handlabeled=length(labels_all);

if sylnum~=sylnum_handlabeled;
    disp('WARNING - hand labeled # of syl does not match actual number, can cause autolabeling frameshift error - check evsonganaly thresholding');
else
disp('GOOD - hand labeled num of syl is correct');
end

lt_jc_fndlbl(vlsorind,vlsorfn,Inds,labels_all,1)


%% NOTES:

% tested onset and offset created by autolabel - they are same as
% evsonganaly.
