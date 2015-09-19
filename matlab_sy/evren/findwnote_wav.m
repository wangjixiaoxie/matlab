function [fvalsstr]=findwnote_wav(batch,NOTE,PRENOTE,POSTNOTE,...
		        TIMESHIFT,FVALBND,NFFT,USEFIT);
%[fvalsstr]=findwnote_wav(batch,NOTE,PRENOTE,POSTNOTE,TIMESHIFT,FVALBND,NFFT,USEFIT);

if (~exist('PRENOTE'))
    PRENOTE='';
elseif (length(PRENOTE)<1)
    PRENOTE='';
end

if (~exist('NFFT'))
    NFFT=1024;
elseif (length(NFFT)<1)
    NFFT=1024;
end

if (~exist('USEFIT'))
    USEFIT=1;
elseif (length(USEFIT)<1)
    USEFIT=1;
end


note_cnt=0;avnote_cnt=0;fcnt=0;
ff=load_batchf(batch);
for ifn=1:length(ff)
    fn=ff(ifn).name;
    fnn=[fn,'.not.mat'];
    if (~exist(fnn,'file'))
        continue;
    end
    disp(fn);
    load(fnn)
    labels = lower(labels);
    labels(findstr(labels,'0'))='-';

    [dat,fs]=wavread(fn);
    if (length(dat)==0)
        disp(['hey no data!']);
        continue;
    end
    fcnt=fcnt+1;

    p=findstr(labels,[PRENOTE,NOTE,POSTNOTE])+length(PRENOTE);
    for ii = 1:length(p)
        ton=onsets(p(ii));toff=offsets(p(ii));

        ti1=ceil((TIMESHIFT + ton*1e-3)*fs);
        if (ti1+NFFT-1<=length(dat))
            note_cnt = note_cnt + 1;
            dattmp=dat([ti1:(ti1+NFFT-1)]);
            fdattmp=abs(fft(dattmp.*hamming(length(dattmp))));

            %get the freq vals in Hertz
            fvals=[0:length(fdattmp)/2]*fs/(length(fdattmp));
	    fdattmp=fdattmp(1:end/2);
            mxtmpvec=zeros([1,size(FVALBND,1)]);
            for kk = 1:size(FVALBND,1)
                tmpinds=find((fvals>=FVALBND(kk,1))&(fvals<=FVALBND(kk,2)));

		NPNTS=2;
                [tmp,pf] = max(fdattmp(tmpinds));
		pf = pf + tmpinds(1) - 1;
                if (USEFIT==1)
                        tmpxv=pf + [-NPNTS:NPNTS];
			tmpxv=tmpxv(find((tmpxv>0)&(tmpxv<=length(fvals))));

			mxtmpvec(kk)=fvals(tmpxv)*fdattmp(tmpxv);
			mxtmpvec(kk)=mxtmpvec(kk)./sum(fdattmp(tmpxv));
		else
                	mxtmpvec(kk) = fvals(pf);
                end
            end

            fvalsstr(note_cnt).fn     = fn;
            fvalsstr(note_cnt).mxvals = [1,mxtmpvec];
            fvalsstr(note_cnt).fdat   = fdattmp;
            fvalsstr(note_cnt).datt   = dattmp;
            fvalsstr(note_cnt).ons    = onsets;
            fvalsstr(note_cnt).offs   = offsets;
            fvalsstr(note_cnt).lbl    = labels;
            fvalsstr(note_cnt).ind    = p(ii);
            fvalsstr(note_cnt).ver    = 5;
            fvalsstr(note_cnt).NPNTS  = NPNTS;
        else
            disp('hey');
        end
    end
end
return;
