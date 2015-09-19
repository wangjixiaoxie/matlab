function [fdatav,freqs,fvalsstr]=get_fft_chunk(batch,NOTE,PRENOTE,POSTNOTE,...
			     	      TIMESHIFT,NFFT,CHANSPEC,ADDX);
%[fdatav,freqs,fvalsstr]=get_fft_chunk(batch,NOTE,PRENOTE,POSTNOTE,...
%			     	      TIMESHIFT,NFFT,CHANSPEC,ADDX);
%[fdat,ff,fv]=get_fft_chunk(bt,NT,PRENT,PSTNT,tbinshft,NFFT,CS);
%


if (~exist('PSTNOTE'))
	PSTNOTE='';
elseif (length(PSTNOTE)<1)
	PSTNOTE='';
end

if (~exist('PRENOTE'))
	PRENOTE='';
elseif (length(PRENOTE)<1)
	PRENOTE='';
end

if (~exist('CHANSPEC'))
	CHANSPEC='obs0';
elseif (length(CHANSPEC)<1)
	CHANSPEC='obs0';
end

if (~exist('ADDX'))
	ADDX=0;
elseif (length(ADDX)==0)
	ADDX=0;
end

ff=load_batchf(batch);
note_cnt=0;avnote_cnt=0;fcnt=0;fdatav=[];
for ifn=1:length(ff)
    fnn=ff(ifn).name;

    if (~exist(fnn,'file'))
        continue;
    end
    disp(fnn);
    load([fnn,'.not.mat']);
    labels = lower(labels);
    labels(findstr(labels,'0'))='-';
    if (ADDX==1)
        ptemp=findstr(fnn,'.cbin');
	if (length(ptemp)==0)
		ptemp=findstr(fnn,'.ebin');
	end
        fnrt=[fnn(1:ptemp(end)-1),'X.rec'];
        rd=readrecf(fnrt);
    else
        rd=readrecf(fnn);
    end

    fcnt=fcnt+1;
    [pthstr,tnm,ext] = fileparts(fnn);
    if (strcmp(ext,'.ebin'))
        [dat,fs]=readevtaf(fnn,CHANSPEC);
    else
        [dat,fs]=evsoundin('',fnn,CHANSPEC);
    end

    if (length(dat)==0)
        disp(['hey no data!']);
        continue;
    end

    p=findstr(labels,[PRENOTE,NOTE,POSTNOTE])+length(PRENOTE);
    for ii = 1:length(p)
        ton=onsets(p(ii));
	toff=offsets(p(ii));

        if (isfield(rd,'ttimes'))
            trigindtmp=find((rd.ttimes>=ton)&(rd.ttimes<=toff));
            if (length(trigindtmp)>0)
                TRIG=1;
                if (isfield(rd,'catch'))
                    ISCATCH=rd.catch(trigindtmp);
                else
                    ISCATCH=-1;
                end
            else
                TRIG=0;
                ISCATCH=-1;
            end
        else
            TRIG=0;
            ISCATCH=-1;
        end

        ti1=ceil((TIMESHIFT + ton*1e-3)*fs);
        if (ti1+NFFT-1<=length(dat))
            note_cnt = note_cnt + 1;
            dattmp=dat([ti1:(ti1+NFFT-1)]);
            fdattmp=abs(fft(dattmp.*hamming(length(dattmp))));


            %get the freq vals in Hertz
            fvals=[0:length(fdattmp)/2]*fs/(length(fdattmp));
            fvals=[fvals,-fvals(end-1:-1:2)];
	    fvals=fvals(1:end/2);fdattmp=fdattmp(1:end/2);

	    if (length(fdatav)==0)
		fdatav=abs(fdattmp);
		freqs=get_fft_freqs(length(fdatav)*2,fs);
	    else
		fdatav=fdatav+abs(fdattmp);
	    end

            [hr,dy,mn,yr]=fn2date(fnn);
            fvalsstr(note_cnt).fn    = fnn;
            fvalsstr(note_cnt).dn    = fn2datenum(fnn);
            fvalsstr(note_cnt).hour  = hr;
            fvalsstr(note_cnt).day   = dy;
            fvalsstr(note_cnt).month = mn;
            fvalsstr(note_cnt).TRIG = TRIG;
            fvalsstr(note_cnt).CATCH = ISCATCH;
            fvalsstr(note_cnt).fdat=fdattmp;
            fvalsstr(note_cnt).datt=dattmp;
            fvalsstr(note_cnt).ons=onsets;
            fvalsstr(note_cnt).offs=offsets;
            fvalsstr(note_cnt).lbl=labels;
            fvalsstr(note_cnt).ind=p(ii);
        else
            disp('hey');
        end
    end
end
fdatav=fdatav./note_cnt;
return;
