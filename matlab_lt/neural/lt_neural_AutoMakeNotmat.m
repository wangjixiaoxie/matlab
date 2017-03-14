function lt_neural_AutoMakeNotmat(batchf)

%% lt 3/9/17 - makes .not mat for all songs in batchf

% need to:
% - give batchf 
% - at least one song has to have segmentation params already (i.e. notmat
% file exists). 
% - all notmat files that currently exist must have same segmentation
% stats.


%%

% 1) get notmat params by looking at all notmat files for this batch
fid=fopen(batchf,'r');

mindurAll = [];
minintAll = [];
thresholdAll = [];
sm_winAll = [];


while (1)
    fn=fgetl(fid);
    
    if (~ischar(fn));
        break;
    end;
    
    if (~exist(fn,'file'))
        continue;
    end
    
    % - does a notmat already exist?
    if exist([fn,'.not.mat'],'file');
        disp(['ALREADY EXISTS: ' fn '.not.mat'])
        TMP = load([fn,'.not.mat']);
        
        mindurAll = [mindurAll TMP.min_dur];
        minintAll = [minintAll TMP.min_int];
        thresholdAll = [thresholdAll TMP.threshold];
        sm_winAll = [sm_winAll TMP.sm_win];
        
    end    
end

assert(~isempty(mindurAll), 'PROBLEM - NO notmat FILE MADE YET. CREATE ONE USING EVSONGANALY');

% ---- COLLECT PARAMS TO USE
MINDUR = unique(mindurAll); assert(length(MINDUR)==1, 'sdafasdf');
MININT = unique(minintAll); assert(length(MININT)==1, 'sdafasdf');
THRESH = unique(thresholdAll); assert(length(THRESH)==1, 'sdafasdf');
SMWIN = unique(sm_winAll); assert(length(SMWIN)==1, 'sdafasdf');
 
% ------------- CREATE ALL NOTMAT FILES THAT DON'T EXIST
fid=fopen(batchf,'r');

while (1)
    fn=fgetl(fid);
    
    if (~ischar(fn));
        break;
    end;
    
    if (~exist(fn,'file'))
        continue;
    end
    
    if (~exist([fn,'.not.mat'],'file'))
        % then create new notmat
        chanspec = '0';
        [dat,Fs,DOFILT,ext]=ReadDataFile(fn,chanspec);
%         [datFull,fs]=evsoundin('',fn,CHANSPEC);
        sm=SmoothData(dat,Fs,1,'hanningfirff');
        sm(1)=0.0; sm(end)=0.0;
        [ons,offs]=SegmentNotesJC(sm,Fs,MININT,MINDUR,THRESH);
        onsets=ons*1e3;offsets=offs*1e3;
        labels = char(ones([1,length(onsets)])*fix('-'));
        
        % == Save notmat and update counter
        fname=fn;
        threshold = THRESH;
        min_dur = MINDUR;
        min_int = MININT;
        sm_win = SMWIN;
        cmd = ['save ',fn,'.not.mat fname Fs labels min_dur min_int ',...
            'offsets onsets sm_win threshold'];
        eval(cmd);
        
        disp(['MADE NEW .not.mat:  ' fn '.not.mat'])
        
        clear sm
    end
    
end
        

        
