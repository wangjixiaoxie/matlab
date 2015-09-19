function [avsp,t,f,avsm]=get_avn(bt,NT,PRETIME,POSTTIME,PRENT,POSTNT,CS,PLTIT);
%[avsp,t,f,avsm]=get_avn(bt,NT,PRETIME,POSTTIME,PRENT,POSTNT,CS,PLTIT);
% bt=batch file
% NT target note
% PRETIME  in seconds
% POSTTIME in seconds
% PRENT  = pre note
% POSTNT = post note
%CS chan spec

if (~exist('PLTIT'))
	PLTIT=0;
elseif (length(PLTIT)==0)
	PLTIT=0;
end

if (~exist('CS'))
	CS='obs0';
elseif (length(CS)==0)
	CS='obs0';
end

if (~exist('PRENT'))
	PRENT='';
elseif (length(PRENT)==0)
	PRENT='';
end

if (~exist('POSTNT'))
	POSTNT='';
elseif (length(POSTNT)==0)
	POSTNT='';
end

if (PLTIT==1)
	figure;
end

fid=fopen(bt,'r');
files=[];cnt=0;
while (1)
	fn=fgetl(fid);
	if (~ischar(fn))
		break;
	end
	cnt=cnt+1;
	files(cnt).fn=fn;
end
fclose(fid);

ontimes=[];

spcnt=0;
for ii=1:length(files)
	fn=files(ii).fn;
	if (exist([fn,'.not.mat'],'file'))
		load([fn,'.not.mat']);
	else
		labels=[];
	end
	labels(findstr(labels,'0'))='-';

	if (ii==1)
		[dat,fs]=evsoundin('',fn,CS);
		NPRE=ceil(PRETIME*fs);
		NPOST=ceil(POSTTIME*fs);
		dat_tmp=zeros([NPRE+NPOST+1,1]);
	end

	disp(fn);
	srchstr=[PRENT,NT,POSTNT];
	pp=findstr(labels,srchstr)+length(PRENT);
	if (length(pp)>0)
		[dat,fs]=evsoundin('',fn,CS);
		
		for jj=1:length(pp)
			tmpon=onsets(pp(jj));
			tmppp=find(abs(tmpon-onsets)<=max([PRETIME,POSTTIME])*1e3);
			ontimes=[ontimes;(onsets(tmppp)-tmpon)*1e-3];
			onind=fix(round(onsets(pp(jj))*1e-3*fs));
			st=onind-NPRE;
			en=onind+NPOST;
			if (st<1)
				stind=abs(st)+2;
				st=1;
			else
				stind=1;
			end
			
			if (en>length(dat))
				enind=length(dat_tmp)-(en-length(dat));
				en=length(dat);
			else
				enind=length(dat_tmp);
			end

			dat_tmp=0*dat_tmp;
			dat_tmp(stind:enind)=dat(st:en);
			
			[sm,sp,t,f]=evsmooth(dat_tmp,fs,10,512,0.8,2,100);
			if(spcnt==0)
				avsp=abs(sp);
				avsm=sm;
				spcnt=1;
			else
				avsp=avsp+abs(sp);
				avsm=avsm+sm;
				spcnt=spcnt+1;
			end
		end
	end
end
t=t-PRETIME;
avsp=avsp./spcnt;
avsm=avsm./spcnt;
disp(num2str(spcnt));
if (PLTIT==1)
	imagesc(t,f,log(abs(avsp)));syn;ylim([0,1e4]);
	hold on;
	plot(ontimes,500+0*ontimes,'k.');
end
return;
