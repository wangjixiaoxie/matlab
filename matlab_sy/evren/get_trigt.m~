function trigs=get_trigt(bt,cntrng,refrac,NFFT,ADDX,TEMPFILE,DO_OR);
%get_trigt(bt,cntrng,refrac,NFFT,ADDX,TEMPFILE,DO_OR);
%
% bt=batch file
% cntrng = struct with MIN, MAX, TH, NOT
% refrac in seconds
% NFFT is the number of points in the templ file
%ADDX if ==1 add the X to the tmp file name to use
% termperary tmp files defualt == 1
% TEMPFILE make a .###X.rec file to avoid overwrite
%    if this is == 1 (defualt)
%DO_OR if ==1 then will use ORs between templates
%
%get_trigt('batch.train',cntrng,.085,128,0,1);
%


if (~exist('ADDX'))
    ADDX=1;
end
if (~exist('TEMPFILE'))
    TEMPFILE=1;
end
if (~exist('NFFT'))
    NFFT=128;
end

if (~exist('refrac'))
    refrac=100e-3;
end

if (~exist('DO_OR'))
    DO_OR=0;
else
    if (length(DO_OR)==0)
        DO_OR=0;
    end
end

for ii=1:length(cntrng)
	if (~isfield(cntrng(ii),'NOT'))
		cntrng(ii).NOT=0;
	end
	if (~isfield(cntrng(ii),'MODE'))
		cntrng(ii).MODE=1;
	end
end

files=load_batchf(bt);

NCNT=length(cntrng);
cnt=zeros([1,NCNT]);
trigs=[];
for ijk=1:length(files)
    fn=files(ijk).name;
    if (~exist(fn,'file'))
        continue;
    end

    pp=findstr(fn,'.cbin');
    if (length(pp)==0)
        pp=findstr(fn,'.ebin');
    end

    if (ADDX==1)
        fnt=[fn(1:pp(end)-1),'X.tmp'];
    else
        fnt=[fn(1:pp(end)),'tmp'];
    end

    if (~exist(fnt,'file'))
        disp('hey');
        continue;
    end
    disp(fnt);

    rdat=readrecf(fn);
    tmpdat=load(fnt);
    tmpdat2=zeros([fix(length(tmpdat)/NCNT),NCNT]);
    for ii=1:NCNT
        tmpdat2(:,ii)=tmpdat(ii:NCNT:end);
    end
    tmpdat=tmpdat2;
    clear tmpdat2;

    fs=rdat.adfreq;
    %TH=1.5;fs=32000;
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
                	cnt(kk)=0;
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

        if (DO_OR==0)
            trig=1;
        else
            trig=0;
        end
        for kk=1:NCNT
            if ((cnt(kk)>=cntrng(kk).MIN)&(cnt(kk)<=cntrng(kk).MAX))
			ntrig=1;
	    else
			ntrig=0;
	    end
	    if (cntrng(kk).NOT==1)
			ntrig=~ntrig;
	    end
            if (DO_OR==0)
                    trig = trig & ntrig;
            else
		    trig = trig | ntrig;
            end
        end

        if (trig)
            if (abs(ii-lasttrig)>refracsam)
                tt=[tt;((ii*NFFT*2/fs)+rdat.tbefore)*1e3];
                lasttrig=ii;
            end
        end
    end
    trigs(ijk).fn=fn;
    trigs(ijk).cntvals=cntvals;

    rdat.ttimes=tt;
    tmp=zeros([length(cntrng),1]);
    for ii=1:length(tmp)
        tmp(ii)=cntrng(ii).TH;
    end
    rdat.thresh=tmp;
    if (TEMPFILE==1)
	fntemp=fn;
	pp=findstr(fn,'.cbin');
	if (length(pp)==0)
		pp=findstr(fn,'.bbin');
	end
	if (length(pp)==0)
		pp=findstr(fn,'.ebin');
	end
	if (length(pp)>0)
		fntemp=[fn(1:pp(end)-1),'X',fn(pp(end):end)];
    		wrtrecf(fntemp,rdat);
	end
     else
    		wrtrecf(fn,rdat);
     end
end
return;
