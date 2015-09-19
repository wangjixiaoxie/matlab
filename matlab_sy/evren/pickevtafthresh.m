function [vals,trigs]=pickevtafthresh(batchf,templ,cntrng,NT,USETEMPF,CHANSPEC,refrac,NFFT);
%[vals,trigs]=pickevtafthresh(batchf,templ,cntrng,NT,USETEMPF,CHANSPEC,refrac,NFFT);
%
% USETEMPF==1 uses tmp file
% refrac -- rerac period in seconds
% NFFT # of points ffted
% 

if (~exist('CHANSPEC'))
    CHANSPEC='obs0';
end

if (~exist('NFFT'))
    NFFT=128;
end

if (~exist('USETEMPF'))
    USETEMPF=0;
end

if (~exist('refrac'))
    refrac = 150.0e-3; % in seconds
end

NCNT=length(cntrng);
cnt=zeros([1,NCNT]);

vals=[];
trigcnt=0;

fid=fopen(batchf,'r');
while (1)
    fn=fgetl(fid);
    if (~ischar(fn))
        break;
    end
    
    fnm=[fn,'.not.mat'];
    if ((~exist(fn,'file'))|(~exist(fnm,'file')))
        continue;
    end
    load(fnm);
    
    if (USETEMPF==1)
        pp=findstr(fn,'.cbin');
        fnt=[fn(1:pp(end)),'tmp'];
        fnr=[fn(1:pp(end)),'rec'];
        if ((~exist(fnt,'file'))|(~exist(fnr,'file')))
            continue;
        end
        tmpdat=load(fnt);
        rd=readrecf(fn);
        fs=rd.adfreq;
        
        tmpdat2=zeros([fix(length(tmpdat)/length(cntrng)),length(cntrng)]);
        for ii = 1:length(cntrng)
            tmpdat2(:,ii)=tmpdat(ii:length(cntrng):end);
        end
        tmpdat=tmpdat2;
        clear tmpdat2;
    else
        [dat,fs] = evsoundin('',fn,CHANSPEC);
        tmpdat   = evtafsim(dat,fs,templ,0,0);
        rd=readrecf(fn);
	tmpdat_t = [1:length(tmpdat)]*2*NFFT/fs;
	postmp=find(tmpdat_t>=rd.tbefore);
	tmpdat=tmpdat(postmp);
    end
    %disp(fn);

    refracsam=ceil(refrac/(2*NFFT/fs));
    lasttrig=-refrac;tt=[];
    cnt=0*cnt;
    for ii = 1:size(tmpdat,1)
        for kk=1:NCNT
            if (tmpdat(ii,kk)<=cntrng(kk).TH)
                cnt(kk)=cnt(kk)+1;
            else
                cnt(kk)=0;
            end
        end

        trig=1;
        for kk=1:NCNT
            if ((cnt(kk)>=cntrng(kk).MIN)&(cnt(kk)<=cntrng(kk).MAX))
                trig = 1 & trig;
            else
                trig = 0;
            end
        end

        if (trig)
            if (abs(ii-lasttrig)>refracsam)
                tt=[tt;((ii*NFFT*2/fs)+rd.tbefore)*1e3];
                lasttrig=ii;
            end
        end
    end
    pp=findstr(labels,NT);
    tmp_out = zeros([1,3]);
    for ii=1:length(tt)
        ppp=find((onsets<=tt(ii))&(offsets>=tt(ii)));
        if (length(ppp)>0)
            if (strcmp(labels(ppp(1)),NT))
                tmp_out(1)=tmp_out(1)+1;
            end
        end
    end
    tmp_out(2) = length(findstr(labels,NT));
    tmp_out(3) = length(tt);
    vals=[vals;tmp_out];
    
    trigcnt=trigcnt+1;
    trigs(trigcnt).fn=fn;
    trigs(trigcnt).tt=tt;
end
fclose(fid);
return;
