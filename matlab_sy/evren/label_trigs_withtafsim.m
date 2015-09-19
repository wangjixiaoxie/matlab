function label_trigs_withtafsim(batch,NT,CSPEC,threshold,cntrng,refrac,NFFT,OVERWRT,min_int,min_dur,sm_win);
% label_trigs_withtafsim(batch,NT,CSPEC,threshold,cntrng,refrac,
%                                     NFFT,OVERWRT,min_int,min_dur,sm_win);
% label_trigs(batch,NT,CSPEC,threshold,OVERWRT,min_int,min_dur,sm_win);
% every trigger is labeled as NT
% if OVERWRT==0 then it will skip files with .not.mat already present
%

if (~exist('OVERWRT'))
    OVERWRT=1;
end

if (~exist('sm_win'))
    sm_win=2.0;
end

if (~exist('min_int'))
    min_int=5.0;
end

if (~exist('min_dur'))
    min_dur=30.0;
end

if (exist('CSPEC'))
    if (~exist('threshold'))
        disp(['Need to specify seg thresh!']);
        return;
    end
end

if (~exist('NFFT'))
    NFFT=128;
end

if (~exist('refrac'))
    refrac=100e-3;
end

fid=fopen(batch,'r');
if (fid==-1)
	disp(['could not find batch file :',batch]);
	return;
end

while (1)
	fn=fgetl(fid);
	if (~ischar(fn))
		break;
	end

	if (~exist(fn,'file'));
		continue;
	end
    if (~exist([fn,'.not.mat'],'file'))
        [dat,Fs]=evsoundin('',fn,CSPEC);
        sm=evsmooth(dat,Fs,0);
        [ons,offs]=SegmentNotes(sm,Fs,min_int,min_dur,threshold);
        onsets=ons*1e3;offsets=offs*1e3;
        fname=fn;
        labels=char(ones([1,length(onsets)])*45);
    else
        if (OVERWRT==0)
            continue;
        end
        load([fn,'.not.mat']);
    end
    
    disp(fn);
    
    NCNT=length(cntrng);
    cnt=zeros([1,NCNT]);

	rd=readrecf(fn);
    rd.ttimes=[];
    
    pp=findstr(fn,'.cbin');
    fnt=[fn(1:pp(end)),'tmp'];
    if (~exist(fnt,'file'))
        continue;
    end
    disp(fn);

    rdat=readrecf(fn);
    tmpdat=load(fnt);
    tmpdat2=zeros([fix(length(tmpdat)/NCNT),NCNT]);
    for ii=1:NCNT
        tmpdat2(:,ii)=tmpdat(ii:NCNT:end);
    end
    tmpdat=tmpdat2;
    clear tmpdat2;

    fs=rdat.adfreq;
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
                tt=[tt;((ii*NFFT*2/fs)+max([rdat.tbefore,0]))*1e3];
                lasttrig=ii;
            end
        end
    end
    rd.ttimes=tt;

	for ii = 1:length(rd.ttimes)
		pp=find((onsets<=rd.ttimes(ii))&(offsets>=rd.ttimes(ii)));
		if (length(pp)>0)
			labels(pp(end))=NT;
		end
	end

	fname=fn;
	cmd = ['save ',fn,'.not.mat fname Fs labels min_dur min_int ',...
	                      'offsets onsets sm_win threshold'];
	eval(cmd);
    clear sm;
end
fclose(fid);
return;
