%% LT 10/1/15 - To replace handchecked mislabeled syls after autolabeling
% After running autolabel it creates a .wav file that contains all the syls
% you labeled.  Use evsonganaly to change labels for mislabeld syls in that
% .wav file. Then run this function to replace those syls in the .cbin song
% files (will replace syls with whatever you label in the .wav files)
function lt_autolabel_FixHandCheckedSyls(fnames, sylnum, vlsorfn, vlsorind)
    
% INPUTS:
% fnames, sylnum, vlsorfn, vlsorind are all outputs from lt_autolabel_function. I.E:
% [fnames, sylnum, vlsorfn, vlsorind]=lt_autolabel_function(Params);
% fnames - filenames of songs in batch
% sylnum - total num of labele dsyls found
% vlsorfn, vlsorind - cell and vector of fiolenames and note positions, for all notes detected
    

% OLD VERSION IGNORE:
%     X1=load('syllwv1.wav.not.mat');
    %     X2=load('syllwv2.wav.not.mat');
    %     XcatLabels=strcat(X1.labels, X2.labels);
    %     [Inds, NewSyl]=regexp(XcatLabels,'\w','start','match');
    % %     newsyl=strjoin(NewSyl,'\0');
    %
    %     lt_jc_fndlbl(vlsorind, vlsorfn, Inds, NewSyl,1)
    
    
    %% ===
    % collect all labels from the .wav file
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
    
    lt_jc_fndlbl(vlsorind,vlsorfn,Inds,labels_all,1);
    


