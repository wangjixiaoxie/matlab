function pj_autoLabel(batchfile,template,cntrng,syl,refrac,ampThresh,min_dur,min_int,EVMODE,NFFT,sm_win)
%% Find all syllables in batch, recognized by template, label with sylLabel. Create 
% labelled .cbin.not.mat files. User should delete existing .cbin.not.mat files before
% running. Runs all steps for each template entered. 
% The arguments templates,cntrngs,and syls should all be lists of strings
% of equal length.

currentDir = pwd;

if (~exist('EVMODE'))
	EVMODE=0;
elseif (length(EVMODE)==0)
	EVMODE=0;
end
if ~exist('sm_win')
    sm_win = 2.0;
end
if (~exist('NFFT'))
    NFFT=128;
end

CHANSPEC='obs0';

    
s = load(template);
fields = fieldnames(s);
templ = s.(fields{1});


% g = load(cntrng);
% fields = fieldnames(g);
% cntrng = g.(fields{1});
cntrng = load(cntrng);

sylLabel = syl;

%% mk_tempf part

NTEMPL=size(templ,2);
fid=fopen(batchfile,'r');
count = 1;
while (1)
	fn=fgetl(fid);
    rdat=readrecf(fn);
    PREDATATIME = rdat.tbefore;
    disp(fn)
	if (~ischar(fn));break;end;
    if (~exist(fn,'file'))
        continue;
    end
    tic
    [dat,fs]=evsoundin('',fn,CHANSPEC);
    timers(count,1) = toc;
    datFull = dat;
    dat=dat(fix(PREDATATIME*fs):end);
    tic
    vals=evtafsim(dat,fs,templ,EVMODE);
    timers(count,2) = toc;

    for ii=1:length(cntrng)     % get_trig2 part
        if (~isfield(cntrng(ii),'NOT'))
            cntrng(ii).NOT=0;
        end
        if (~isfield(cntrng(ii),'MODE'))
            cntrng(ii).MODE=1;
        end
    end

    NCNT=length(cntrng);
    cnt=zeros([1,NCNT]);
    trigs=[];
    pp=findstr(fn,'.cbin');
    tmpdat=vals;
    fs=rdat.adfreq;
    refracsam=ceil(refrac/(2*NFFT/fs));
    lasttrig=-refrac;tt=[];
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
    rdat.ttimes=tt;
    tic
    if (~exist([currentDir,'/',fn,'.not.mat'],'file'))
        if(size(dat)==0)
            continue
        end
        [sm]=SmoothData(datFull,fs,1,'hanningfirff');
        sm(1)=0.0;sm(end)=0.0;
        [ons,offs]=SegmentNotes(sm,fs,min_int,min_dur,ampThresh);
        onsets=ons*1e3;offsets=offs*1e3;
        labels = char(ones([1,length(onsets)])*fix('-'));
    else
        load([fn,'.not.mat']);
    end
	for ii = 1:length(rdat.ttimes)
		pp=find((onsets<=rdat.ttimes(ii))&(offsets>=rdat.ttimes(ii)));
		if (length(pp)>0)
			labels(pp(end))=sylLabel;
        end
    end
	fname=fn;
    Fs = fs;
	cmd = ['save ',fn,'.not.mat fname Fs labels min_dur min_int ',...
	                      'offsets onsets sm_win ampThresh'];
	eval(cmd);
    clear sm;
    count = count + 1;
    timers(count,4) = toc;
end
end

%fclose(fid);
%return;


