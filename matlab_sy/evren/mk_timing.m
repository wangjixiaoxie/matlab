function [av,tvals,st,cnt]=mk_timing(bt,NT,PRET,PSTT,fs,CS);
% [av,tvals,sm,cnt]=mk_timing(bt,NT,PRET,PSTT,fs,CS);
% bt=batch file
% NT-note ('a')
% PRET-PRETIME in seconds (1)
% PSTT-POSTTIME in seconds (1)
% fs - sample freq in hz (32e3)
%

if (~exist('CS'))
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
st=av;
cnt=av;
for ii=1:length(fns)
	disp([fns(ii).fn]);
	load([fns(ii).fn,'.not.mat']);
    pp=findstr(labels,NT);
    [dat,fs]=evsoundin('',fns(ii).fn,CS);
    %if (~exist([fns(ii).fn,'.sm'],'file'))
        [sm]=evsmooth(dat,fs,0.0,512,0.8,2,1500);
    %else
    %    eval(['load -mat ',fns(ii).fn,'.sm']);
    %end
    for jj=1:length(pp)
        ton=onsets(pp(jj))*1e-3;
		inds=floor((ton-PRET)*fs)+[0:length(av)-1];
		ppp=find((inds>0)&(inds<=length(sm)));

		av(ppp) = av(ppp)  + sm(inds(ppp)).';
		st(ppp) = st(ppp)  +(sm(inds(ppp)).').^2;
		cnt(ppp)= cnt(ppp) + 1;
	end
end
tvals=[1:length(av)]/fs - PRET;
av=av./cnt;
st=sqrt(st./cnt - av.^2);
return;
