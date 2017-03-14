%% LT 2/18/17 - added output with segment of analog data that is extra channels (i.e. not audio)

%% LT 6/15/15 - added catch song output (not just catch trial as before)
%% LT 12/21/14 -
% added:
% 1) ttimes
% 2) iscatch (syl level)
% 3) trigger note (i.e. for evtav4)
% to output

% 4) ADDX - if 1, then uses X.rec. if 0, then uses normal .recs

%% LT 7/24/14 - Modified from findwnote2tw -
% 1) annotations
% 2) allows evtafv4 or evtafamp for getting time from filename (see end)
% 3) allows option to not run fft and just extract sound data.  the fft is largely useless unless you have a timeshift that makes sense. Otherwise use pitchcountour or evtaffreq (at template match site);
% 4) added fvalsstr(note_cnt).NotePos
% 5) Timeshift should be 0.016. that's becuase the first sample will be 16ms after 1st time point, since NFFT of 1024 is 32ms (fs=32000), so midpoint of that window is 16ms. (sliding window, 8000hz)


%% LT 7/7/ removed lower from line 47 so keeps capital labels
function [fvalsstr]=findwnote2tw_v4_LT(batchnotmat,NOTE,PRENOTE,TIMESHIFT,FVALBND,NFFT,...
    ADDNOTMAT,USEFIT,CHANSPEC,calcFFT,ADDX)
% NEW PINTERP VERSION
% [fvalsstr]=findwnote2(batchnotmat,NOTE,PRENOTE,TIMESHIFT,FVALBND,NFFT,...
%                                                       ADDNOTMAT,USEFIT,CHANSPEC);
%
%
%
% LT notes:
% NFFT - number of samples to take for pitch analysis (i.e. 32,000hz)
% calcFFT = 0 (no) or 1 (yes);
% USEFIT - interpolates value of fr at power peak.




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
if (~exist('calcFFT')) % Default is to not run fft.
    calcFFT=0;
end


fid=fopen(batchnotmat,'r');
% AVN= zeros([257,111]);
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
    
    if exist('ADDX');
        rd     = readrecf(fnrt,ADDX);
    else
        rd     = readrecf(fnrt);
    end
    
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

    %% ============== read data from other channels if exist
    dat_otherchans = {};
    if rd.nchan>1
        for i=2:rd.nchan
       [dat_otherchans{i-1},fs]=ReadDataFile(fnn(1:pos(end)-1),num2str(i-1)); 
        end 
    end
    
    %%
    
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
        
        % GETS data, including raw sound data.
        ton=onsets(p(ii));toff=offsets(p(ii)); % only get to this point if pre and note both match.
        
        % Note whether trigger occured at this note.
        NoteNum=[];
        CatchTrial=[];
        if (isfield(rd,'ttimes'))
            if (length(find((rd.ttimes>=ton)&(rd.ttimes<=toff)))>0)
                TRIG=1;
                
                % if trig, then was it catch trial?
                TrigInds=find((rd.ttimes>=ton)&(rd.ttimes<=toff));
                if length(TrigInds)==1;
                    CatchTrial=rd.catch(TrigInds); % 1 if catch trial, 0 if not
                    
                    % get notenum
                    yy=findstr(rd.pbname{TrigInds},'Templ');
                    NoteNum=str2double(rd.pbname{TrigInds}(yy+8));
                    
                elseif length(TrigInds)>1;
                    
                    CatchTrial=rd.catch(TrigInds);
                    
                    yy=findstr(rd.pbname{TrigInds(1)},'Templ');
                    NoteNum=str2double(rd.pbname{TrigInds(1)}(yy+8));
                    disp(['WARNING - multiple triggers for one syl in song ' fname])
                    
                end
                
                
            else
                TRIG=0;
            end
        else
            TRIG=0;
        end
        
        note_cnt = note_cnt + 1;
        pttmpvec=zeros([1,size(FVALBND,1)]);
        
        ti1=ceil((TIMESHIFT + ton*1e-3)*fs);
        dattmp=dat([ti1:(ti1+NFFT-1)]);
        
        % ----------------------------------------------
        % FREQUENCY CALCULATION - spectrum taken at one time bin (sepcified
        % by timeshift)
        if calcFFT==1;
            fdattmp=abs(fft(dattmp)).^2;
            fvals=[0:length(fdattmp)/2]*fs/(length(fdattmp));
            fvals=[fvals,-fvals(end-1:-1:2)];
            mxtmpvac=zeros([size(FVALBND,1),2]);
            tmpinds = find((fvals>=FVALBND(1))&(fvals<=FVALBND(2)));
            
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
                    tmpyv = sqrt(abs(fdattmp(tmpinds(tmpxv))));
                    tmpyv = tmpyv./sum(tmpyv);
                    
                    pf = tmpxv*tmpyv+tmpinds(1)-1;
                else
                    pf = pf + tmpinds(1) - 1;
                end
            else
                pf = pf + tmpinds(1) - 1;
            end
            mxtmpvec = (pf-1)*(fvals(2)-fvals(1)); % multiply bin # by size of bins in frequency
            
        end
        % -----------------------------------------
        
        % OUTPUT STUFF:
        datenum=fn2datenum_eftafv4_lt(fnrt);
        %         [hr,dy,mn,yr]=fn2date(fnrt) % old, but doesn't work for both
        %         evtaf v4 and amp
        fvalsstr(note_cnt).fn    = fnn;
        %         fvalsstr(note_cnt).hour  = hr;
        %         fvalsstr(note_cnt).day   = dy;
        %         fvalsstr(note_cnt).month = mn;
        
        fvalsstr(note_cnt).datenum=datenum;
        fvalsstr(note_cnt).TRIG = TRIG;
        fvalsstr(note_cnt).datt=dattmp; % sound data
        fvalsstr(note_cnt).ons=onsets;
        fvalsstr(note_cnt).offs=offsets;
        fvalsstr(note_cnt).lbl=labels;
        fvalsstr(note_cnt).NotePos=p(ii);
        
        try
            fvalsstr(note_cnt).ttimes=rd.ttimes;
        catch err
        end
        
        % what template notes were used during triggering?
        if length(rd.ttimes)>0
            for jj=1:length(rd.ttimes);
                yy=findstr(rd.pbname{jj},'Templ');
                fvalsstr(note_cnt).NoteNums(jj)=str2double(rd.pbname{jj}(yy+8));
            end
            
        end
        
        % what were the catch statuses of those triggers?
        if length(rd.catch)>0;
            fvalsstr(note_cnt).CatchList=rd.catch;
        end
    
        fvalsstr(note_cnt).CatchTrial=CatchTrial;
        fvalsstr(note_cnt).NoteNum=NoteNum;

        % --- Song was catch song?
        fvalsstr(note_cnt).CatchSong=rd.iscatch;
        
        try % i.e. maybe didn't do FFT
            fvalsstr(note_cnt).mxvals = [1,mxtmpvec];
            fvalsstr(note_cnt).fdat=fdattmp; % frequency data (squared fft).
        catch err
        end
        
        % ==== IF EXTRA ANALOG CHANS, THEN SAVE
        if rd.nchan>1
            dat_otherchans_tmp = {};
            for i=2:rd.nchan
                % -- confirm is correct length
                assert(length(dat_otherchans{i-1}) == length(dat), 'asdfasdf');
                
                % -- EXTRACT CORRECT TIMEPOINT
                dat_otherchans_tmp{i-1} = dat_otherchans{i-1}([ti1:(ti1+NFFT-1)]);
                            end
            % --- save
            fvalsstr(note_cnt).dat_otherchans=dat_otherchans_tmp; % sound data
        end
        
    end
end

fclose(fid);
% if (avnote_cnt>0)
%     AVN=AVN/avnote_cnt;
% end
return;
