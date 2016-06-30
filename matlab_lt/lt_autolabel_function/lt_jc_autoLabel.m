function lt_jc_autoLabel(batchfile,templatefile,cntrng,syl,refrac,ampThresh,min_dur,min_int,overwrite_notmat)
%% LT 4/9/15 - modified (changes pending)
% Things to modify: ALL DONE
% 1) can hit a little outside onset and offset.
% 2) can run multiple times with different templates and accumulate labels
% (not overwriting each time). did not do. instead made option to overwrite
% or append old not mat. effectively same thing
% 3) can choose to label +/- the target. (e.g. string)
% 4) resaves batch file adding AutoLabeled to name
% 5) saves old notmat files in backup folder


% New inputs:
% overwrite_notmat=1 (will overwrite old not mat. i.e. will save copy of
% notmat in backup folder (only for songs in this batch) and create new
% notmat.
% overwrite_notmat=0 (if notmat exists, will add onto it. i.e. keeps old
% labels, unless that label conflicts with new label). also saves a
% backup copy of old notmat
% syl is structure telling how to fill in syls.
% syl.targ='b' label trigger
% syl.pre='cc' string preceding trigger
% syl.post='dd' string post trigger
%     (note, if not detected notes (e.g. edge of song), then will do as much as possible).
        % template is actually filename (e.g. template.dat). it should be
        % saved in dir one up.

%% Old description. Find all syllables in batch, recognized by template, label with sylLabel. Create 
% labelled .cbin.not.mat files. User should delete existing .cbin.not.mat files before
% running. Runs all steps for each template entered. 
% The arguments templates,cntrngs,and syls should all be lists of strings
% of equal length.
% template: use matrix with each column a templ
% cntrng: make structure with:
% 
%     cntrng(1).MIN=minCount;
%     cntrng(1).MAX=maxCount;
%     %true/false logic, true->note=0
%     cntrng(1).NOT=0;
%     %evtafmode=1; birdtafmode=0;
%     cntrng(1).MODE=1;
%     %threshold
%     cntrng(1).TH=FFTthresh;
%     %and/or logic with other templates.
%     cntrng(1).AND=0;
%     cntrng(1).BTMIN=0;

%% INITIATE PARAMS
sm_win = 2; % used in case will save new .not.mat files
NFFT = 128; % frequency bins for fft, means window is 256 samples.
CHANSPEC='obs0'; % for opening labview files
EVMODE=0; % mode for nromalizing templates when performing evtafsim.  Keep it 0 which means use new method.

SylBoundsAdd=16; % in ms, how much time to remove from onset and add to offset in permissive window for labeling? I make it 0.016s (2 fft bins)
% this is useful for the time when a syl is fragmented (e.g. problem
% with thresholding), or trig is outside onset or offset by a little.

% How many syls in string used for labeling?
NumLabelsPre=length(syl.pre);
NumLabelsPost=length(syl.post);
LabelStr=[syl.pre syl.targ syl.post];

% if NumLabelsPost==0 && NumLabelsPre==0;
%     sylLabel = syl.targ;
%     old_way=1; % only one target, slightly quicker. (about 90% time).
% else
%     old_way=0;
% end

%% make folder to save backups of old notmats
tstamp=lt_get_timestamp(0);

OldNotMatdir=['OldNotMat_AL_moved_' tstamp];
mkdir(OldNotMatdir);

%% LOAD TEMPLATE, assume it is in above dir.
currentDir = pwd;
template=load(['../' templatefile]);


%% mk_tempf part

fid=fopen(batchfile,'r');
count = 1;
while (1)
    fn=fgetl(fid);
    disp(fn)
    
    if (~ischar(fn));
        break;
    end;
    
    if (~exist(fn,'file'))
        continue;
    end
    
    rdat=readrecf(fn);
    PREDATATIME = rdat.tbefore;
    
    [dat,fs]=evsoundin('',fn,CHANSPEC);
    
    if(size(dat)==0)
        continue
    end
    
    datFull = dat;
    dat=dat(fix(PREDATATIME*fs):end);
    vals=evtafsim(dat,fs,template,EVMODE); % these are .tmp values
    
    
    for ii=1:length(cntrng)     % get_trig2 part
        if (~isfield(cntrng(ii),'NOT'))
            cntrng(ii).NOT=0;
        end
        if (~isfield(cntrng(ii),'MODE'))
            cntrng(ii).MODE=1; % default mode is not evtaf
        end
    end

    NCNT=length(cntrng);
    cnt=zeros([1,NCNT]);
    trigs=[];
    pp=findstr(fn,'.cbin');
    tmpdat=vals;
    refracsam=ceil(refrac/(2*NFFT/fs));
    
    % initiate
    lasttrig=-refrac;
    tt=[];
    cnt=0*cnt;
    for kk=1:NCNT
        if (cntrng(kk).MODE==0)
            cnt(kk)=cntrng(kk).MAX+1;
        end
    end
	
    cntvals=zeros([NCNT,size(tmpdat,1)]);
    for ii = 1:size(tmpdat,1)
        for kk=1:NCNT
            if (tmpdat(ii,kk)<=cntrng(kk).TH)
                if (cntrng(kk).MODE==1)
                    cnt(kk)=cnt(kk)+1;
                else
                    if (cnt(kk)>=cntrng(kk).BTMIN)
                        cnt(kk)=0;
                    else
                        cnt(kk)=cnt(kk)+1;
                    end
                end
            else
                if (cntrng(kk).MODE==0)
                    cnt(kk)=cnt(kk)+1;
                else
                    cnt(kk)=0;
                end
            end
        end
	cntvals(:,ii)=cnt.';
        for kk=1:NCNT
            if ((cnt(kk)>=cntrng(kk).MIN)&(cnt(kk)<=cntrng(kk).MAX))
                ntrig=1;
            else
                ntrig=0;
            end
            if (cntrng(kk).NOT==1)
                ntrig=~ntrig;
            end
            if (kk==1)
                trig=ntrig;
            else
                if (cntrng(kk-1).AND==1)
                    trig = trig & ntrig;
                else
                    trig = trig | ntrig;
                end
            end
        end
        
        if (trig)
            if (abs(ii-lasttrig)>refracsam)
                tt=[tt;((ii*NFFT*2/fs)+rdat.tbefore)*1e3];
                lasttrig=ii;
            end
        end
    end
    
    % output (trigger times):
    rdat.ttimes=tt;
    
    % does a notmat already exist?
    if (~exist([currentDir,'/',fn,'.not.mat'],'file')); % does not already exist
        
        % then create new notmat
        [sm]=SmoothData(datFull,fs,1,'hanningfirff');
        sm(1)=0.0;sm(end)=0.0;
        [ons,offs]=SegmentNotes(sm,fs,min_int,min_dur,ampThresh);
        onsets=ons*1e3;offsets=offs*1e3;
        labels = char(ones([1,length(onsets)])*fix('-'));
    else
        % first save backup of old notmat
        eval(['!cp ' fn '.not.mat ' OldNotMatdir]);
        
        % second, perform one of two options
        if overwrite_notmat==1;
            % delete old notmat
            eval(['!rm ' fn '.not.mat'])
            
            % create new notmat
            [sm]=SmoothData(datFull,fs,1,'hanningfirff');
            sm(1)=0.0;sm(end)=0.0;
            [ons,offs]=SegmentNotes(sm,fs,min_int,min_dur,ampThresh);
            onsets=ons*1e3;offsets=offs*1e3;
            labels = char(ones([1,length(onsets)])*fix('-'));
            
        else
            % load old notmat -
            load([fn,'.not.mat']);
        end
    end
    
    
    % Mark new labels
    for ii = 1:length(rdat.ttimes)
        pp=find(((onsets-SylBoundsAdd)<=rdat.ttimes(ii))&((offsets+SylBoundsAdd)>=rdat.ttimes(ii)));
        %         find(((onsets)<=rdat.ttimes(ii))&((offsets)>=rdat.ttimes(ii)))
        if (length(pp)>0);
            %                 labels(pp(1))=sylLabel; % greater than 1 is almost certainly multiple trigs, unless SylBoundsAdd was set larger than smallest gap dur, take 1st syl with trig.
            try % try to stick flanking sequence in. if fails, then fits as much as possible on either side.
                labels(pp(1)-NumLabelsPre:pp(1)+NumLabelsPost)=LabelStr; % greater than 1 is almost certainly multiple trigs, unless SylBoundsAdd was set larger than smallest gap dur, take 1st syl with trig.
            catch err
                disp(['error: song ' num2str(ii) ' could not fit all label string, will fill with what I can']);
                
                % find out how much of label can stick in there.
                if pp(1)-1<NumLabelsPre;
                    Npre=pp(1)-1;
                else
                    Npre=NumLabelsPre;
                end
                if length(onsets)-pp(1)<NumLabelsPost;
                    Npost=length(onsets)-pp(1);
                else
                    Npost=NumLabelsPost;
                end
                
                % fill in labels as can
                labels(pp(1)-Npre:pp(1)+Npost)=LabelStr(NumLabelsPre-Npre+1:end-(NumLabelsPost-Npost));
            end
        end
    end
    
    % Save notmat
    fname=fn;
    Fs = fs;
    threshold = ampThresh;
    cmd = ['save ',fn,'.not.mat fname Fs labels min_dur min_int ',...
        'offsets onsets sm_win threshold'];
	eval(cmd);
    clear sm;
    count = count + 1;
    
end

% Save new version of batch file with "autolabeled" in name
batch_new=[batchfile '.AutoLabeled'];
eval(['!cp ' batchfile ' ' batch_new])

disp('Done autolabeling!')

end