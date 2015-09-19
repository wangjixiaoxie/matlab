function [avsp,t,f]=get_avn(bt,NT,PRETIME,POSTTIME,PRENT,POSTNT,CS);
%[avsp,t,f]=get_avn(bt,NT,PRETIME,POSTTIME,PRENT,POSTNT,CS);
% bt=batch file
% NT target note
% PRETIME  in seconds
% POSTTIME in seconds
% PRENT  = pre note
% POSTNT = post note
%CS chan spec

if (~exist('CS'))
	CS='obs0';
end
if (~exist('PRENT'))
	PRENT='';
end
if (~exist('POSTNT'))
	POSTNT='';
end


fid=fopen(bt,'r');
files=[];cnt=0;
while (1)
	fn=fgetl(fid);
	if (~ischar(fn));
		break;
	end
	cnt=cnt+1;
	files(cnt).fn=fn;
end
fclose(fid);

spcnt=0;
for ii=1:length(files);
	fn=files(ii).fn;
	if (exist([fn,'.not.mat'],'file'));
		load([fn,'.not.mat']);
	else
		labels=[];
	end


	if (ii==1);
		[dat,fs]=evsoundin('',fn,CS);
		NPRE=ceil(PRETIME*fs);
		NPOST=ceil(POSTTIME*fs);
		dat_tmp=zeros([NPRE+NPOST+1,1]);
	end

	disp(fn);
	srchstr=[PRENT,NT,POSTNT];
	pp=findstr(labels,srchstr)+length(PRENT);
	if (length(pp)>0);
		[dat,fs]=evsoundin('',fn,CS);
		
		for jj=1:length(pp);
			onind=fix(round(onsets(pp(jj))*1e-3*fs));
			st=onind-NPRE;
			en=onind+NPOST;
			if (st<1);
				stind=abs(st)+2;
				st=1;
			else
				stind=1;
			end
			
			if (en>length(dat));
				enind=length(dat_tmp)-(en-length(dat));
				en=length(dat);
			else
				enind=length(dat_tmp);
			end

			dat_tmp=0*dat_tmp;
			dat_tmp(stind:enind)=dat(st:en);
			
			[sm,sp,t,f]=evsmooth(dat_tmp,fs,0,512,0.8,2,100);
			if(spcnt==0)
				avsp=abs(sp);
				spcnt=1;
			else
				avsp=avsp+abs(sp);
				spcnt=spcnt+1;
			end
		end
	end
end
t=t-PRETIME;
avsp=avsp./spcnt;
disp(num2str(spcnt));
return;
