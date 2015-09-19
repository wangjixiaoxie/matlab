function trigs=get_trigt(bt,cntrng,refrac,NFFT,ADDX,TEMPFILE);
%trigs=get_trigt2(bt,cntrng,refrac,NFFT,ADDX,TEMPFILE);
%
% bt=batch file
% cntrng = struct with MIN, MAX, TH, NOT, BTMIN, AND (==1 means do and)
% refrac in seconds
% NFFT is the number of points in the templ file
%ADDX if ==1 add the X to the tmp file name to use
% termperary tmp files defualt == 1
% TEMPFILE make a .###X.rec file to avoid overwrite
%    if this is == 1 (defualt)
%DO_OR if ==1 then will use ORs between templates
%
% ADDED THE NEW COUNTER VERSION WERE IS THE COUNTER HAS
% TO BE >=N and val<TH for the count to continue goin up


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

%if (~exist('DO_OR'))
%    DO_OR=0;
%else
%    if (length(DO_OR)==0)
%        DO_OR=0;
%    end
%end

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
        disp('hey, X.tmp does not exist - should run mk_tempf to make those');
        continue;
    end
    disp(fnt);
    
    rdat=readrecf(fn);
    tmpdat=load(fnt); % loads temp file
    tmpdat2=zeros([fix(length(tmpdat)/NCNT),NCNT]);
    for ii=1:NCNT % converts temp file to multi-column format
        placeholder=tmpdat(ii:NCNT:end);
        tmpdat2(:,ii)=placeholder(1:size(tmpdat2,1));
    end
    tmpdat=tmpdat2;
    clear tmpdat2;
    
    fs=rdat.adfreq;
    %TH=1.5;fs=32000;
    refracsam=ceil(refrac/(2*NFFT/fs));  % how many 0.008s samples are there in refract
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
        
        %if (DO_OR==0)
        %    trig=1;
        %else
        %    trig=0;
        %end
        
        % have to separate into bracket groups here
        
        for kk=1:NCNT
            if ((cnt(kk)>=cntrng(kk).MIN)&(cnt(kk)<=cntrng(kk).MAX))
                ntrig=1;
            else
                ntrig=0;
            end
            if (cntrng(kk).NOT==1)
                ntrig=~ntrig;
            end
            %if (DO_OR==0)
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
        % now have all bracket groups, what is the logic between those groups?
        
        % Finally, single output: trig or no trig?
        if (trig) % output time is the time right after the end of the trigger bin.
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
