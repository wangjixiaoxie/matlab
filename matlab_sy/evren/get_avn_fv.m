function [avsp,t,f,avsm]=get_avn_fv(fv,NT,PRENT,POSTNT,PLTIT);
% [avsp,t,f,avsm]=get_avn_fv(fv,PRETIME,POSTTIME,PRENT,POSTNT,PLTIT);

if (~exist('PLTIT'))
	PLTIT=0;
elseif (length(PLTIT)==0)
	PLTIT=0;
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

ontimes=[];
spcnt=0;
for ii=1:length(fv)
	labels=fv(ii).lbl;
	dat=fv(ii).datt;
	ind=fv(ii).ind;

	inds=ind+[-length(PRENT):length(POSTNT)];
	srchstr=[PRENT,NT,POSTNT];
	if (strcmp(labels(inds),srchstr))
		[sm,sp,t,f]=evsmooth(dat,32e3,0);
		if (spcnt==0)
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
avsp=avsp./spcnt;
avsm=avsm./spcnt;
disp(num2str(spcnt));
if (PLTIT==1)
	imagesc(t,f,log(abs(avsp)));syn;ylim([0,1e4]);
	hold on;
	plot(ontimes,500+0*ontimes,'k.');
end
return;
