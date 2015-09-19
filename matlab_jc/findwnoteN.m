%findwnote6
%modified to add STIMTRIG
%modified to add STIMISCATCH

%modified to deal with new evtaf4 file format.


function [fvalsstr]=findwnote9(batch,NOTE,PRENOTE,POSTNOTE,...
		        TIMESHIFT,FVALBND,NFFT,USEFIT,CHANSPEC,ADDX,ONSETWINDOW, SOUNDSTR, STIMSTR);
% NEW PINTERP VERSION
% [fvalsstr]=findwnote2(batch,NOTE,PRENOTE,POSTNOTE,...
%		TIMESHIFT,FVALBND,NFFT,USEFIT,CHANSPEC,ADDX);
%
%fv=findwnote4(bt,NT,PRENT,PSTNT,tbinshft,fbins,NFFT,1,'obs0');

% soundstr='wav'
% stimstr='rig'

fvalsstr=[];
if (~exist('PRENOTE'))
    PRENOTE='';
elseif (length(PRENOTE)<1)
    PRENOTE='';
end

if(~exist('ONSETADJ'))
    ONSETADJ=0;
end
    

if (~exist('CHANSPEC'))
    CHANSPEC='obs0';
elseif (length(CHANSPEC)<1)
    CHANSPEC='obs0';
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

if (~exist('ADDX'))
    ADDX=0;
elseif (length(ADDX)<1)
    ADDX=0;
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
    load(fnn);
% %     if(exist('mod_onset'))
% %         onsets=mod_onset;
% %     end
    
    labels = lower(labels);
    labels(findstr(labels,'0'))='-';
    if (ADDX==1)
        ptemp=findstr(fnn,'.cbin');
	if (length(ptemp)==0)
        	ptemp=findstr(fnn,'.cbin');
	end
		
        fnrt=[fnn(1:ptemp(end)-1),'X.rec'];
    else
	fnrt=fn;
    end
    rd = readrecf(fnrt);

    [pthstr,tnm,ext] = fileparts(fn);
    if (strcmp(ext,'.ebin'))
        [dat,fs]=readevtaf(fn,CHANSPEC);
    else
        [dat,fs]=evsoundin('',fn,CHANSPEC);
    end

    if (length(dat)==0)
        disp(['hey no data!']);
        continue;
    end
    fcnt=fcnt+1;

    p=findstr(labels,[PRENOTE,NOTE,POSTNOTE])+length(PRENOTE);
    for ii = 1:length(p)
        if(length(onsets)==length(labels))
            ton=onsets(p(ii))-ONSETWINDOW(1);toff=onsets(p(ii))-ONSETWINDOW(2);
        
            if (isfield(rd,'ttimes'))
            trigindtmp=find((rd.ttimes>=ton)&(rd.ttimes<=toff));
                if (length(trigindtmp)>0)
                    trigon=0;
                    soundtrigon=0;
                    for trigind=1:length(trigindtmp)
                        if(~isempty(findstr(STIMSTR,rd.pbname{trigindtmp(trigind)})))
                            trigon=1;
                            trigvl=trigind
                        end
                        if(~isempty(findstr(SOUNDSTR,rd.pbname{trigindtmp(trigind)})))
                            soundtrigon=1;
                            soundtrigvl=trigind
                        end

                    end
                    if(soundtrigon)
                        SOUNDTIME=rd.ttimes(trigindtmp(soundtrigvl))-onsets(p(ii));
                    else
                        SOUNDTIME=[];
                    end

                    if(trigon)
                        STIMTRIG=1;
                        STIMTIME=rd.ttimes(trigindtmp(trigvl))-onsets(p(ii));
                        SOUNDTRIG=0;
                        if(isfield(rd,'catch'))
                        	STIMCATCH=rd.catch(trigindtmp(trigvl));
                            SOUNDCATCH=0;
                        
                        else
                            SOUNDTRIG=1;
                            STIMTRIG=0;
                            SOUNDCATCH=rd.catch(trigindtmp(trigvl));
                            STIMCATCH=0;
                        end
                    
                
                     else
                        TRIG=0;
                        STIMTRIG=0
                        SOUNDTRIG=0;
                        ISCATCH=-1;
                        SOUNDCATCH=-1
                        STIMCATCH=-1
                     end
                else
                    TRIG=0;
                    STIMTRIG=0;
                    SOUNDTRIG=0;
                    SOUNDTIME=[];
                    ISCATCH=-1
                    SOUNDCATCH=-1
                    STIMCATCH=-1
                end
                end
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

                NPNTS=10;
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

            [hr,dy,mn,yr]=fn2date(fn);
            fvalsstr(note_cnt).fn     = fn;
            fvalsstr(note_cnt).hour   = hr;
            fvalsstr(note_cnt).day    = dy;
            fvalsstr(note_cnt).month  = mn;
            fvalsstr(note_cnt).mxvals = [1,mxtmpvec];
            fvalsstr(note_cnt).STIMTRIG   = STIMTRIG;
            fvalsstr(note_cnt).STIMCATCH  = STIMCATCH;
            fvalsstr(note_cnt).SOUNDTRIG   = SOUNDTRIG;
            fvalsstr(note_cnt).SOUNDCATCH  = SOUNDCATCH;
            fvalsstr(note_cnt).STIMTIME=STIMTIME;
            fvalsstr(note_cnt).SOUNDTIME=SOUNDTIME;
            
            fvalsstr(note_cnt).fdat   = fdattmp;
            fvalsstr(note_cnt).datt   = dattmp;
            fvalsstr(note_cnt).ons    = onsets;
            fvalsstr(note_cnt).offs   = offsets;
            fvalsstr(note_cnt).lbl    = labels;
            fvalsstr(note_cnt).ind    = p(ii);
            fvalsstr(note_cnt).ver    = 4;
            fvalsstr(note_cnt).NPNTS  = NPNTS;
        else
            disp('hey');
            end
        else
            fid2=fopen('~/output.txt','a+');
                fprintf(fid2, '%s\n',fn);
                fclose(fid2);
        end
    end
end
return;
