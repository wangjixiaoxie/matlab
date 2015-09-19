function [AVN,fvals,hvals,mxvals,vecs]=findwnote(batchnotmat, ...
    NOTE,PRENOTE,POSTNOTE,TBINSHIFT,FBINBND,ADDNOTMAT,USEFIT,CHANSPEC,NFFT);
% NEW PINTERP VERSION
%[AVN,fvals,hvals,mxvals,vecs]=findwnote(batchnotmat,NOTE,PRENOTE,POSTNOTE,...
%                       TBINSHIFT,FBINBND,ADDNOTMAT,USEFIT,CHANSPEC,NFFT);
%

if (~exist('PRENOTE'))
    PRENOTE='';
elseif (length(PRENOTE)<1)
    PRENOTE='';
end
if (~exist('CHANSPEC'))
    CHANSPEC='obs0r';
end
if (~exist('NFFT'))
    NFFT=512;
end
if (~exist('ADDNOTMAT'))
    ADDNOTMAT=0;
end
if (~exist('USEFIT'))
    USEFIT=0;
end

fid=fopen(batchnotmat,'r');
vecs=[];mxvals=[];
AVN= zeros([257,111]);
note_cnt=0;
avnote_cnt=0;
fcnt=0;
while (1)
    fnn=fgetl(fid);
    
    if (~ischar(fnn))
        break;
    end
    
    
    if (ADDNOTMAT)
        fnn = [fnn,'.not.mat'];
    end
    
    if (~exist(fnn,'file'))
        continue;
    end
    
    disp(fnn);load(fnn);
    pos = findstr(fnn,'.not.mat');
    fnrt = fnn(1:pos(end)-1);
    labels=lower(labels);
    rd=readrecf(fnrt);

    fcnt=fcnt+1;
    [pthstr,tnm,ext] = fileparts(fnn(1:pos(end)-1));
    if (strcmp(ext,'.ebin'))
        [dat,fs]=readevtaf(fnn(1:pos(end)-1),CHANSPEC);
    else
        [dat,fs]=evsoundin('',fnn(1:pos(end)-1),CHANSPEC);
    end
    %if (fix(fs)~=32e3)
    %    tt=([1:length(dat)]-1)/fs;
    %    dat=interp1(tt,dat,[0:1.0./32e3:tt(end)].','linear');
    %    fs=32e3;
    %end
    if (length(dat)==0)
        continue;
    end
    [sm,sp,t,f]=evsmooth(dat,fs,0);

    %p=findstr(labels,NOTE);
    p=findstr(labels,[PRENOTE,NOTE,POSTNOTE])+length(PRENOTE);
    for ii = 1:length(p)
        %if (length(PRENOTE)>0)
        %    if (p(ii)>length(PRENOTE))
        %        lblinds=[(p(ii)-length(PRENOTE)):(p(ii)-1)];
        %        if (~strcmp(PRENOTE,labels(lblinds)))
        %            continue;
        %        end
        %    end
        %end
       
        ton=onsets(p(ii));toff=offsets(p(ii));
        [tmp,ind1] = min(abs(ton*1e-3  - t));
        [tmp,ind2] = min(abs(toff*1e-3 - t));
        
        if (isfield(rd,'ttimes'))
            if (length(find((rd.ttimes>=ton)&(rd.ttimes<=toff)))>0)
                TRIG=1;
            else
                TRIG=0;
            end
        else
            TRIG=0;
        end

        note_cnt = note_cnt + 1;
        mnfreqs=zeros([length(TBINSHIFT),size(FBINBND,1)]);
        mxvals=[];
        ptvals=zeros([length(TBINSHIFT),size(FBINBND,1)]);
        for ijk = 1:length(TBINSHIFT)
            mxtmpvec=[];
            
            pttmpvec=zeros([1,size(FBINBND,1)]);
            delt=t(2)-t(1);
            ti1=floor((TBINSHIFT(ijk)*delt + ton*1e-3)*fs)+1;
            dattmp=dat([ti1:(ti1+NFFT-1)]);
            fdattmp=abs(fft(dattmp)).^2;
            fdattmp=xcorr(fdattmp,'coeff');
            fdattmp = fdattmp(length(dattmp):end);

            for kk = 1:size(FBINBND,1)
                tmpinds = [FBINBND(kk,1):FBINBND(kk,2)];
                tmpsp=abs(sp(tmpinds,ind1+TBINSHIFT(ijk)-1));
                MXV  = max(max(tmpsp));
                [pf,pt] = find(tmpsp==MXV);
                if (USEFIT==1)
                    if ((pf>2)&(pf<size(tmpsp,1)-1))
                        pkfit=polyfit(pf+[-2:2],tmpsp(pf+[-2:2]).',2);
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

                fdattmp2=fdattmp(tmpinds);
                [y,i]=max(fdattmp2);
                if ((i<3)|(i>length(fdattmp2)-2))
                    xmax = i;
                else
                    %xmax=pinterp(i+[-2:2].',fdattmp2(i+[-2:2].'));
                    polyvals=polyfit(i+[-2:2].',fdattmp2(i+[-2:2].'),2);
                    xmax=-polyvals(2)./2./polyvals(1);
                end
                pttmpvec(kk) = xmax + tmpinds(1) - 1;
            end
            mxvals = [mxvals;ijk,mxtmpvec];
            ptvals(ijk,:) = pttmpvec;
        end
	if (strcmp(CHANSPEC,'w'))
		hr=0;
		dy=0;
		mn=0;
	else
		[hr,dy,mn,yr]=fn2date(fnrt);
	end
        fvals(note_cnt).fn    = fnn;
        fvals(note_cnt).hour  = hr;
        fvals(note_cnt).day   = dy;
        fvals(note_cnt).month = mn;
        fvals(note_cnt).mnfreqs = mnfreqs;
        fvals(note_cnt).mxvals  = mxvals;
        fvals(note_cnt).TRIG = TRIG;
        fvals(note_cnt).pterp=ptvals;
        fvals(note_cnt).fdat=fdattmp;
        fvals(note_cnt).datt=dattmp;
        fvals(note_cnt).ons=onsets;
        fvals(note_cnt).lbl=labels;

        tmpvec = abs(sp(:,ind1+TBINSHIFT-1));
        for ll=1:size(tmpvec,2)
            tmpvec(:,ll) = tmpvec(:,ll)-min(tmpvec(:,ll));
            tmpvec(:,ll) = tmpvec(:,ll)./max(tmpvec(:,ll));
        end
        vecs = [vecs,tmpvec];

        if (((ind1+100)<=size(sp,2))&(ind1>=11))
            AVN= AVN+ abs(sp(:,ind1+[-10:100]));
            avnote_cnt = avnote_cnt+1;
        end
    end
end
fclose(fid);
if (avnote_cnt>0)
    AVN=AVN/avnote_cnt;
end

% do the histograms for all the TBINSHIFTS
hvals = zeros([size(AVN,1),length(TBINSHIFT)+1]);
hvals(:,1) = [1:size(AVN,1)].';
if (length(mxvals)>0)
    for ii = 1:length(TBINSHIFT)
        pp = find(mxvals(:,1)==ii);
        tmpdat = mxvals(pp,2:end);
        tmpdat = reshape(tmpdat,[size(tmpdat,1)*size(tmpdat,2),1]);
        [b,a] = hist(tmpdat,hvals(:,1));
        hvals(:,ii+1) = b;
    end
else
    disp(['no mxvals?']);
end
return;
