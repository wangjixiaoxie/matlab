%% LT 7/7/ removed lower from line 47 so keeps capital labels
function [fvalsstr]=findwnote2tw(batchnotmat,NOTE,PRENOTE,TIMESHIFT,FVALBND,NFFT,...
                                                       ADDNOTMAT,USEFIT,CHANSPEC);
% NEW PINTERP VERSION
% [fvalsstr]=findwnote2(batchnotmat,NOTE,PRENOTE,TIMESHIFT,FVALBND,NFFT,...
%                                                       ADDNOTMAT,USEFIT,CHANSPEC);
%
%
%

if (~exist('PRENOTE'))
    PRENOTE='';
elseif (length(PRENOTE)<1)
    PRENOTE='';
end
if (~exist('CHANSPEC'))
    CHANSPEC='obs0';
end
if (~exist('NFFT'))
    NFFT=1024;
end
if (~exist('ADDNOTMAT'))
    ADDNOTMAT=0;
end
if (~exist('USEFIT'))
    USEFIT=0;
end

fid=fopen(batchnotmat,'r');
AVN= zeros([257,111]);
note_cnt=0;avnote_cnt=0;fcnt=0;
while (1)
    fnn=fgetl(fid);
    if (~ischar(fnn))
        break;
    end
    if (ADDNOTMAT)
        fnn = [fnn,'.not.mat'];
    end
    if (~exist(fnn,'file'))
        disp('ERROR, FINDNOTE: FILE DOES NOT EXIST')
        continue;
    end
    disp(fnn);load(fnn)
    pos    = findstr(fnn,'.not.mat');
    fnrt   = fnn([1:pos(end)-1]);
%     labels = lower(labels);
    rd     = readrecf(fnrt);

    fcnt=fcnt+1;
    [pthstr,tnm,ext] = fileparts(fnn(1:pos(end)-1));
    if (strcmp(ext,'.ebin'))
        [dat,fs]=readevtaf(fnn(1:pos(end)-1),CHANSPEC);
    else
        [dat,fs]=evsoundin('',fnn(1:pos(end)-1),CHANSPEC);
    end

    if (length(dat)==0)
        disp(['hey no data!']);
        continue;
    end
    %[sm,sp,t,f]=evsmooth(dat,fs,0);
    
    % ------------------------
    % Pick out syls that have both pre and note match.
    p=findstr(labels,NOTE);
    for ii = 1:length(p)
        if (length(PRENOTE)>0)
            if (p(ii)>length(PRENOTE)) % i.e. if p is not first labeled syl.
                lblinds=[(p(ii)-length(PRENOTE)):(p(ii)-1)];
                if (~strcmp(PRENOTE,labels(lblinds))) % if note before syl is NOT what you want (prenote).
                    continue;
                end
            end
        end
    % ------------------------

        ton=onsets(p(ii));toff=offsets(p(ii)); % only get to this point if pre and note both match.

        % Note whether trigger occured at this note.
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
        pttmpvec=zeros([1,size(FVALBND,1)]);

        %Timeshift is if you don't want to look at start of note?%
        ti1=ceil((TIMESHIFT + ton*1e-3)*fs);
        dattmp=dat([ti1:(ti1+NFFT-1)]);
        fdattmp=abs(fft(dattmp)).^2;
        
        %get the freq vals in Hertz
        fvals=[0:length(fdattmp)/2]*fs/(length(fdattmp));
        fvals=[fvals,-fvals(end-1:-1:2)];
        mxtmpvac=zeros([size(FVALBND,1),2]);
        for kk = 1:size(FVALBND,1)
            tmpinds = find((fvals>=FVALBND(kk,1))&(fvals<=FVALBND(kk,2)));

            [tmp,pf] = max(fdattmp(tmpinds));
            if (USEFIT==1)
                if ((pf>=3)&(pf<=length(tmpinds)-2))
                    %pkfit=polyfit(pf+[-2:2],fdattmp(pf+[-2:2]).',2);
                    %pkfit=-pkfit(2)./2./pkfit(1);
                    %[pkfit,tmp]=pinterp(pf+[-2:2].',log(fdattmp(pf+[-2:2])));
                    
                    %pkfit=polyfit(pf+[-2:2].',log(fdattmp(pf+[-2:2])),2);
                    %pkfit=roots(pkfit(1:2).*[2:-1:1]);
                    
                    %pf = pkfit + tmpinds(1) - 1;
                    
                    tmpxv = pf+[-2:2];
                    tmpyv = sqrt(abs(fdattmp(tmpinds(pf+[-2:2]))));
                    tmpyv = tmpyv./sum(tmpyv);

                    pf = tmpxv*tmpyv+tmpinds(1)-1;
                else
                    pf = pf + tmpinds(1) - 1;
                end
            else
                pf = pf + tmpinds(1) - 1;
            end
            mxtmpvec(kk) = (pf-1)*(fvals(2)-fvals(1));
        end
	
	[hr,dy,mn,yr]=fn2date(fnrt);
        fvalsstr(note_cnt).fn    = fnn;
        fvalsstr(note_cnt).hour  = hr;
        fvalsstr(note_cnt).day   = dy;
        fvalsstr(note_cnt).month = mn;
        fvalsstr(note_cnt).mxvals = [1,mxtmpvec];
        fvalsstr(note_cnt).TRIG = TRIG;
        fvalsstr(note_cnt).fdat=fdattmp;
        fvalsstr(note_cnt).datt=dattmp;
        fvalsstr(note_cnt).ons=onsets;
        fvalsstr(note_cnt).offs=offsets;
        fvalsstr(note_cnt).lbl=labels;
    end
end
fclose(fid);
if (avnote_cnt>0)
    AVN=AVN/avnote_cnt;
end
return;
