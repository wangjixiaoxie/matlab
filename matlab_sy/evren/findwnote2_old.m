function [AVNOTE,mxhi,vecs,tod]=findwnote2(batchnotmat,NOTE,...
                                                 PRENOTE,ADDNOTMAT);
% VERSION FOR o9bk98
%[AVNOTE,mxhi,vecs]=findwnote2(batchnotmat,NOTE,PRENOTE,ADDNOTMAT);
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

fid=fopen(batchnotmat,'r');
vecs=[];
mxlow=[];mxhi=[];
AVNOTE = zeros([257,111]);
cnt1=0;cnt2=0;
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
    labels=lower(labels);
	pos = findstr(fn,'.not.mat');

	[dat,fs]=readevtaf(fn(1:pos(end)-1),'0r');
	[sm,sp,t,f]=evsmooth(dat,fs,0.01);

	p=findstr(labels,NOTE);
	fbinbnd = fix([20,50;
	          51,90;
		  91,130;
		  131,150]);
		  
	tbinshft=[25,30,35,39];
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
		
		for ijk = 1:length(tbinshft)
			mxtmpvec=[];
			for kk = 1:size(fbinbnd,1)
				tmpinds = [fbinbnd(kk,1):fbinbnd(kk,2)];
				tmpsp=abs(sp(tmpinds,ind1+tbinshft(ijk)-1));
				MXV   = max(max(tmpsp));
				[pf,pt] = find(tmpsp==MXV);
				pf = pf + fbinbnd(kk,1) - 1;
				mxtmpvec = [mxtmpvec,pf];
			end
			mxhi = [mxhi;ijk,mxtmpvec];
		end

		tmpvec = sum(abs(sp(:,ind1+tbinshft-1).')).';
		tmpvec = tmpvec-min(tmpvec);
		tmpvec = tmpvec./max(tmpvec);
		vecs = [vecs,tmpvec];

        if ((ind1+100)<=size(sp,2))
            AVNOTE = AVNOTE + abs(sp(:,ind1+[-10:100]));
            cnt1 = cnt1+1;
        end
	end
end
fclose(fid);
if (cnt1>0)
	AVNOTE=AVNOTE/cnt1;
end
return;
