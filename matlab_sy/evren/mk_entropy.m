function [av,tvals,st,cnt]=mk_entropy(bt,NT,PRET,PSTT,fs,CS);
% [av,tvals,st,cnt]=mk_entropy(bt,NT,PRET,PSTT,fs,CS);
% bt=batch file
% NT-note ('a')
% PRET-PRETIME in seconds (1)
% PSTT-POSTTIME in seconds (1)
% fs - sample freq in hz (32e3)
% CS - channel spec (default 'obs0'

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
for ii=1:length(fns)
	disp([fns(ii).fn]);
	load([fns(ii).fn,'.not.mat']);
 	pp=findstr(labels,NT);
	if (length(pp)==0)
		continue;
	end
    	[dat,fst]=evsoundin('',fns(ii).fn,CS);
        [sm,sp,t,f]=evsmooth(dat,fst,0.0);
	ppp=find((f>=500)&(f<=1e4));
	sp=abs(sp(ppp,:));
	sp=sp./(ones([size(sp,1),1])*sum(sp));
	ent=sum(-sp.*log2(sp));
	ent_t=t;
	
	if (ii==1)
		dt = t(2)-t(1);
		fs = 1.0./dt;
		NPRE=fix(PRET*fs);
		NPST=fix(PSTT*fs);
		av = zeros([1,(NPST+NPRE+1)]);
		st = zeros([1,(NPST+NPRE+1)]);
		cnt= zeros([1,(NPST+NPRE+1)]);
	end

    	for jj=1:length(pp)
        	ton=onsets(pp(jj))*1e-3;
		[y,i]=min(abs(ton-ent_t));
		
		inds=i+[-NPRE:NPST];
		ppp=find((inds>0)&(inds<=length(sm)));

		av(ppp)  =  av(ppp) + ent(inds(ppp));
		st(ppp)  =  st(ppp) + (ent(inds(ppp)).^2);
		cnt(ppp) = cnt(ppp) + 1;
	end
end
tvals=[-NPRE:NPST]/fs;
av=av./cnt;
st=sqrt((st./cnt )-(av.^2));
return;
