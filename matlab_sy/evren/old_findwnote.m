function [AVNOTE,histvals,mxvals,freqvals,vecs]=findwnote(batchnotmat,NOTE,PRENOTE,...
                        TBINSHIFT,FBINBND,ADDNOTMAT,USEFIT);
% GENERAL VERSION
%[AVNOTE,histvals,mxvals,freqvals,vecs]=findwnote(batchnotmat,NOTE,PRENOTE,...
%                       TBINSHIFT,FBINBND,ADDNOTMAT,USEFIT);
%
%DIND=4;
%NOTE='b';
%batchnotmat='batchnotmat';
%ADDNOTMAT=0;

if (~exist('PRENOTE'))
	PRENOTE='';
elseif (length(PRENOTE)<1)
	PRENOTE='';
end
if (~exist('ADDNOTMAT'))
	ADDNOTMAT=0;
end
if (~exist('USEFIT'))
    USEFIT=0;
end
ADDNOTMAT,USEFIT
fid=fopen(batchnotmat,'r');
vecs=[];mxvals=[];
AVNOTE = zeros([257,111]);
note_cnt=0;
avnote_cnt=0;
fcnt=0;
while (1)
	fn=fgetl(fid);
	if (~ischar(fn))
		break;
	end
	if (ADDNOTMAT)
		fn = [fn,'.not.mat'];
	end
	if (~exist(fn,'file'))
		continue;
	end
	disp(fn);load(fn);
	pos = findstr(fn,'.not.mat');
	fnrt = fn(1:pos(end)-1);
	labels=lower(labels);
    rd=readrecf(fnrt);

	fcnt=fcnt+1;
    [pthstr,tnm,ext] = fileparts(fn(1:pos(end)-1));
    if (strcmp(ext,'.ebin'))
        [dat,fs]=readevtaf(fn(1:pos(end)-1),'0r');
    else
        [dat,fs]=evsoundin('',fn(1:pos(end)-1),'obs0r');
    end
    if (length(dat)==0)
        continue;
    end
	[sm,sp,t,f]=evsmooth(dat,fs,0.01);

	p=findstr(labels,NOTE);
    for ii = 1:length(p)
        if (length(PRENOTE)>0)
            if (p(ii)>length(PRENOTE))
                lblinds=[(p(ii)-length(PRENOTE)):(p(ii)-1)];
                if (~strcmp(PRENOTE,labels(lblinds)))
                    continue;
                end
            end
        end
        ton=onsets(p(ii));toff=offsets(p(ii));
		[tmp,ind1] = min(abs(ton*1e-3  - t));
		[tmp,ind2] = min(abs(toff*1e-3 - t));
      
        if (length(find((rd.ttimes>=ton)&(rd.ttimes<=toff)))>0)
            TRIG=1;
        else
            TRIG=0;
        end
		
		note_cnt = note_cnt + 1;
		mnfreqs=zeros([length(TBINSHIFT),size(FBINBND,1)]);
        mxvals=[];
		for ijk = 1:length(TBINSHIFT)
			mxtmpvec=[];
			for kk = 1:size(FBINBND,1)
				tmpinds = [FBINBND(kk,1):FBINBND(kk,2)];
				tmpsp=abs(sp(tmpinds,ind1+TBINSHIFT(ijk)-1));
				MXV  = max(max(tmpsp));
				[pf,pt] = find(tmpsp==MXV);
                if (USEFIT==1)
                    if ((pf>1)&(pf<size(tmpsp,1)))
                        pkfit=polyfit(pf+[-1:1],tmpsp(pf+[-1:1]),2);
                        pkfit=-pkfit(2)./2./pkfit(1);
                        pf = pkfit + FBINBND(kk,1) - 1;
                    else
                        pf = pf + FBINBND(kk,1) - 1;
                    end
                else
                    pf = pf + FBINBND(kk,1) - 1;
                end
                mxtmpvec = [mxtmpvec,pf];
				mnfreqs(ijk,kk) = sum(tmpsp.'.*tmpinds)./sum(tmpsp);
			end
			mxvals = [mxvals;ijk,mxtmpvec];
		end
		[hr,dy,mn,yr]=fn2date(fnrt);
		freqvals(note_cnt).hour  = hr;
		freqvals(note_cnt).day   = dy;
		freqvals(note_cnt).month = mn;
		freqvals(note_cnt).mnfreqs = mnfreqs;
        freqvals(note_cnt).mxvals  = mxvals;
        freqvals(note_cnt).TRIG = TRIG;

		tmpvec = abs(sp(:,ind1+TBINSHIFT-1));
		for ll=1:size(tmpvec,2)
			tmpvec(:,ll) = tmpvec(:,ll)-min(tmpvec(:,ll));
			tmpvec(:,ll) = tmpvec(:,ll)./max(tmpvec(:,ll));
		end
		vecs = [vecs,tmpvec];

		if ((ind1+100)<=size(sp,2))
			AVNOTE = AVNOTE + abs(sp(:,ind1+[-10:100]));
			avnote_cnt = avnote_cnt+1;
		end
	end
end
fclose(fid);
if (avnote_cnt>0)
	AVNOTE=AVNOTE/avnote_cnt;
end

% do the histograms for all the TBINSHIFTS
histvals = zeros([size(AVNOTE,1),length(TBINSHIFT)+1]);
histvals(:,1) = [1:size(AVNOTE,1)].';
for ii = 1:length(TBINSHIFT)
	pp = find(mxvals(:,1)==ii);
	tmpdat = mxvals(pp,2:end);
	tmpdat = reshape(tmpdat,[size(tmpdat,1)*size(tmpdat,2),1]);
	[b,a] = hist(tmpdat,histvals(:,1));
	histvals(:,ii+1) = b;
end
return;
