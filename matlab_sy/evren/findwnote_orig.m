function [AVNOTE,mxvals,vecs]=findwnote(batchnotmat,NOTE,FBINBNDS,TBINBNDS,...
                                                 PRENOTE,ADDNOTMAT);
%[AVNOTE,mxhi,vecs]=
%      findwnote(batchnotmat,NOTE,FBINBNDS,TBINBNDS,PRENOTE,ADDNOTMAT);
%
%DIND=4;
%NOTE='b';
%batchnotmat='batchnotmat';
%FBINMIN=72;
%ADDNOTMAT=0;

if (~exist('PRENOTE'))
	PRENOTE='';
elseif (length(PRENOTE)<1)
	PRENOTE='';
end
if (~exist('ADDNOTMAT'))
	ADDNOTMAT=0;
end

fid=fopen(batchnotmat,'r');
vecs=[];mxvals=[];
AVNOTE = zeros([257,111]);
note_cnt=0;
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

	[dat,fs]=readevtaf(fn(1:pos(end)-1),'0r');
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

		tmpsp = abs(sp([FBINMIN:end],[ind1:ind2]));
		MXV = max(max(tmpsp));
		[pf,pt] = find(tmpsp==MXV);

		% lines below are for finding the averages
		%AVAL = abs(sp(FBINMIN:85,pt).');
		%AVAL = AVAL*f(FBINMIN:85)./sum(AVAL);
		%mxhi = [mxhi;pt AVAL MXV];
		
		pf = pf + FBINMIN - 1;
		mxhi = [mxhi;pt pf MXV];

		%tmpvec = abs(sp(:,ind1+pt(1)-1));
		tmpvec = sum(abs(sp(:,ind1+pt(1)-2:ind1+pt(1)-1).'));
		tmpvec = tmpvec-min(tmpvec);
		tmpvec = tmpvec./max(tmpvec);
		vecs = [vecs,tmpvec.'];

		tmpsp = sp([1:FBINMIN],[ind1:ind2]);
		MXV = max(max(tmpsp));
		[pf,pt] = find(tmpsp==MXV);
		mxlow = [mxlow;pt f(pf) MXV];

		AVNOTE = AVNOTE + abs(sp(:,ind1+[-10:100]));
		note_cnt = note_cnt+1;
	end
end
fclose(fid);
if (note_cnt>0)
	AVNOTE=AVNOTE/note_cnt;
end
return;
