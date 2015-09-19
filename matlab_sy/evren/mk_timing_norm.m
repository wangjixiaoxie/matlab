function [av,tvals,st,cnt]=mk_timing_norm(bt,NT,PRET,PSTT,fs,CS,PRENT,PSTNT);
% [av,tvals,st,cnt]=mk_timing_norm(bt,NT,PRET,PSTT,fs,CS,PRENT,PSTNT);
% bt=batch file
% NT-note ('a')
% PRET-PRETIME in seconds (1)
% PSTT-POSTTIME in seconds (1)
% fs - sample freq in hz (32e3)
% CS - channel spec (default 'obs0'
% PRENT - pre context
% PSTNT - post context

if (exist('PRENT','var')==0)
	PRENT='';
elseif (length(PRENT)==0)
	PRENT='';
end
if (exist('PSTNT','var')==0)
	PSTNT='';
elseif (length(PSTNT)==0)
	PSTNT='';
end

if (exist('CS','var')==0)
	CS='obs0';
elseif (length(CS)==0)
	CS='obs0';
end

if (~exist('NT'))
	NT='a';
elseif (length(NT)==0)
	NT='a';
end
if (~exist('PRET'))
	PRET=1;
elseif (length(PRET)==0)
	PRET=1;
end
if (~exist('PSTT'))
	PSTT=1;
elseif (length(PSTT)==0)
	PSTT=1;
end
if (~exist('fs'))
	fs=32e3;
elseif (length(fs)==0)
	fs=32e3;
end


fid=fopen(bt,'r');
fns=[];
while (1)
	fn=fgetl(fid);
	if (~ischar(fn));
		break
	end
	fns(length(fns)+1).fn=fn;
end
fclose(fid);

%PRET=1.0;PSTT=1.0;fs=32e3;
av=zeros([1,fix((PSTT+PRET)*fs)]);
st=zeros([1,fix((PSTT+PRET)*fs)]);
cnt=zeros([1,fix((PSTT+PRET)*fs)]);
for ii=1:length(fns)
	disp([fns(ii).fn]);
	load([fns(ii).fn,'.not.mat']);
 	pp=findstr(labels,[PRENT,NT,PSTNT])+length(PRENT);
	if (length(pp)==0)
		continue;
	end
	%if (length(pp)==1)
	%	disp(['hey']);
	%end
    	[dat,fs]=evsoundin('',fns(ii).fn,CS);
    	%if (~exist([fns(ii).fn,'.sm'],'file'))
        [sm]=evsmooth(dat,fs,0.0);
    	%else
    	%    eval(['load -mat ',fns(ii).fn,'.sm']);
    	%end

	tst=floor(onsets(pp(1))*1e-3*fs);
	ten=ceil(offsets(pp(end))*1e-3*fs);

	tst=floor(onsets(1)*1e-3*fs);
	ten=ceil(offsets(end)*1e-3*fs);
	normv=mean(sqrt(sm(tst:ten)));

    	for jj=1:length(pp)
        	ton=onsets(pp(jj))*1e-3;
		inds=floor((ton-PRET)*fs)+[0:length(av)-1];
		ppp=find((inds>0)&(inds<=length(sm)));

		av(ppp)  =  av(ppp) + sqrt(sm(inds(ppp)).')./normv;
		st(ppp)  =  st(ppp) + (sm(inds(ppp)).')./normv./normv;
		cnt(ppp) = cnt(ppp) + 1;
	end
end
tvals=[1:length(av)]/fs - PRET;
av=av./cnt;
st=sqrt((st./cnt )-(av.^2));
return;
