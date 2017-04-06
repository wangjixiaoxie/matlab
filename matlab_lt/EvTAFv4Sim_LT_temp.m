function [AllSongsAllTrigsCollected]=EvTAFv4Sim_LT(batchf,configfile,ChanSpec)
%% LT 3/15/17 - handles manual trigs
% previously issue was looked for tmeplate in rec file, but did not find
% since it was a manual trig. now puts a "nan" instead of template number.
% see line ~215

%% LT 6/5/15 - Polished. UP TO DATE SUMMARY:
% Does 3 things:
% 1) Gets offline detects (for all note numbers) using evtaf sim.  Saves those data to a) rec file (X.rec) and to output structure
% 2) Get online info (from online created .rec file) and saves to output structure.
% 3) DOES NOT save that information in new X.rec file.  It makes the rd structure, but I have not modified writerecf function to save the additional stuff.


%% LT 4/22/15 - now collects all data (i.e. all detects, all songs) in an output structure
% AllSongsAllTrigsCollected(n), where n are songs, and each song has data.
% checked if slows timing, and only went from 72 to 73.5 seconds for a set
% of songs.


%% LT 11/25/14 -
% Changes:
% 1) around line 113, previously was giving errors, since if any song file did not have triggers, it would have an empty matrix, leading to error. Changed by having if clause./
% 2) function spits out ND structure
% 3) EvTAFv4_LT instead of EvTAFv4
% 4) Instead of writing actual triggers (i.e. pass Freq and ampl thresh) to
% rec file, will write all triggers (i.e. as long as template matches, and
% repeats, delay and so forth match) to rec file.  This is in line with
% other versions of this analysis (i.e. get_trigt2_v4_LT).
% NOTE: saves the real trigs to another field in the rec file, and also
% saves Frequency values for real and all trigs.

% NOTE - DOES NOT SAVE CATCH TRIALS - see wrtrecf_LT_evtafv4 notes for
% reason

%%
%ND=EvTAFv4Sim(batchf,configfile,ChanSpec);
% Loads every cbin file in the batchfile -> batchf and simulates
%  EvTAFv4 triggers.  All EvTAFv4 params are loaded from the file
% listed in configifle.  This must be new binary
% config file format.  Note there is an old version of the config file
% format that will not work here but this should only effect me (evren)
% and tim warren.  ChanSpec is the channel you want to use as the
% sound/song data.  usualy format 'obs0' or 'obs2r' etc.
%
%output - write X.rec files with the simulated trigger times in them
% so if the cbin filename is : bk57w35_170809_093631.355.cbin
% this program will output a file called : bk57w35_170809_093631.355X.rec
% this way any original rec files do not get written over.
% you can read this recfile by using readrecf which returns a strucutre
% with all the file info.  Just set the second input ADDX to 1 and it will
% add the X in X.rec for you automatically.  example:
% rd=readrecf('bk57w35_170809_093631.355.cbin',1);
% then rd.ttimes is the trigger times in msec
% also useful rd.trignote is which of the notes triggered meaning which
% index of Note Detection Info brought about the trigger
%
%ND is an array of strucutures.  it contains all kinds of
%information about the evtaf triggering. It corresponds to the Note
%Detection Info strucutre in the main EvTAFv4 LabVIEW window
% if you want to see the trigger times

[ND,OP]=ReadEvTAFv4ConfigFile(configfile);

% do some checks to make sure everything is ok
for ii=1:length(ND)
    if isempty(ND(ii).Templ)
        disp(['Templates are empty for note detect info index = ',num2str(ii)]);
        if exist(ND(ii).TemplFile,'file')
            disp(['Load the templates from the template file : ',ND(ii).TemplFile]);
            ND(ii).Templ=load(ND(ii).TemplFile);
        else
            disp(['Could not file the Template File indicated : ',ND(ii).TemplFile]);
            disp(['Probably a path mismatch']);
            disp(['Please choose the right template file']);
            [FileName,PathName,FilterIndex] = uigetfile('*',['Template for Note Index ',num2str(ii)]);
            ND(ii).TemplFile=[FileName,PathName];
            ND(ii).Templ=load(ND(ii).TemplFile);
        end
    end
    
    %verify template normalizations
    for jj=1:size(ND(ii).Templ,2)
        ND(ii).Templ(1:6,jj)=0;
        ND(ii).Templ(:,jj)=ND(ii).Templ(:,jj)./max(ND(ii).Templ(:,jj));
    end
    
    % get the number of data points that have to be read in based on the
    % size of the templates
    if (ii==1)
        NFFT=2*size(ND(ii).Templ,1);
    else
        if (NFFT~=2*size(ND(ii).Templ,1))
            disp(['ERROR - EvTAF Limitation all templates have the be the same length in freq bins']);
            return;
        end
    end
    
    if (length(ND(ii).CntRng)~=size(ND(ii).Templ,2))
        disp(['Counter Range values and Templates are of different sizes for Note index : ',num2str(ii)]);
        disp(['Cannot continue']);
        return;
    end
end


inputfiles=[];
fid=fopen(batchf,'r');
while (1)
    fn=fgetl(fid);
    if (ischar(fn))
        if (exist(fn,'file'))
            inputfiles(length(inputfiles)+1).fn=fn;
        else
            disp(['Could not fine file : ',fn]);
            disp('Skipping it');
        end
    else
        break;
    end
end
fclose(fid);



AllSongsAllTrigsCollected=struct; 

for IFile=1:length(inputfiles)
    %loop over every file in the batchfile
    fn=inputfiles(IFile).fn;
    disp(fn);
    [dat,fs]=evsoundin('',fn,ChanSpec); % load up the raw data
    
    PREBUFIND=ceil(OP.FileBufferLeng*fs);
    
    %skip the data in the prebuffer - Cannot trigger off this data since
    %it is just held in the buffer and no sound playback happens if the box
    %is 'quiet'
    
    [ND_EvSimDone,TemplMatchVals]=EvTAFv4_LT(ND,OP,dat,fs);
    
    SimTrigInfoActual=[];
    SimTrigInfoAll=[];
    
    FreqValsActual=[];
    FreqValsAll=[];
    
    
    % ==== SAVE INFORMATION ABOUT THIS SONG
    
    AllSongsAllTrigsCollected(IFile).filename=fn;
    AllSongsAllTrigsCollected(IFile).datenum_filestart_SecResol=fn2datenum_eftafv4_lt(fn); % seconds resol, when file start, not including pre-buffer
    AllSongsAllTrigsCollected(IFile).data_OfflineAllDetects=[];
    AllSongsAllTrigsCollected(IFile).header={'trigtime_ms','notenum','FF','Amp','TRIG_offline'};
    
    % -------------------
    
    % LT added saving frequency values.
    for iND=1:length(ND_EvSimDone)
        
        % LT - old version, which saves only real triggers to rec file.
        % THESE ARE ACTUALLY DEPENDENT ON FF THRESHOLD OF CONFIG FILE!! NOT
        % NECESSARILY REAL.  ALSO, THESE DO NOT CARE ABOUT NOTE GROUPS
        % (Functionality of Evtaf_v4_LT).  Accurate online trigs are in rec
        % file.  This here might not be totally accurate for those 2
        % reasons above.  i.e. this is not used:
        SimTrigInfoActual=[SimTrigInfoActual;...
            ND_EvSimDone(iND).TriggerTimes*1e3,ones(size(ND_EvSimDone(iND).TriggerTimes))*(iND-1), ND_EvSimDone(iND).FFVals]; % col 1 are trigger times, col 2 are note IDs (0,...)
        
        
        
        % LT - new version, saves all Triggers to rec file
        SimTrigInfoAll=[SimTrigInfoAll;...
            ND_EvSimDone(iND).AllTrigs.TriggerTimes*1e3, ones(size(ND_EvSimDone(iND).AllTrigs.TriggerTimes))*(iND-1), ND_EvSimDone(iND).AllTrigs.FFVals]; % col 1 are trigger times, col 2 are note IDs (0,...)
        
        
        % --------------------------
        % Collect data for all detects across all songs.

        AllSongsAllTrigsCollected(IFile).data_OfflineAllDetects=[AllSongsAllTrigsCollected(IFile).data_OfflineAllDetects; ...
            ND_EvSimDone(iND).AllTrigs.TriggerTimes*1e3, ones(size(ND_EvSimDone(iND).AllTrigs.TriggerTimes))*(iND-1), ND_EvSimDone(iND).AllTrigs.FFVals, ND_EvSimDone(iND).AllTrigs.AmpVals];

        % ----------------------------
            
    end
    
    % old version, LT changed to below
    %         [sortv,sorti]=sort(SimTrigInfo(:,1));
    %         SimTrigInfo=SimTrigInfo(sorti,:);
    %         rd=readrecf(fn);
    %         rd.ttimes=SimTrigInfo(:,1);
    %         rd.trignote=SimTrigInfo(:,2);
    
    % Modified again, see below (to incorporate both actual and All trigs)
    %     if ~isempty(SimTrigInfo); % i.e. this song contains at least one trig.
    %         [sortv,sorti]=sort(SimTrigInfo(:,1));
    %         SimTrigInfo=SimTrigInfo(sorti,:);
    %         rd=readrecf(fn);
    %         rd.ttimes=SimTrigInfo(:,1);
    %         rd.trignote=SimTrigInfo(:,2);
    %     else
    %         rd=readrecf(fn);
    %         rd.ttimes=[];
    %         rd.trignote=[];
    %     end
    
    
    rd=readrecf(fn);
    
    % ==== SAVE INFORMATION ABOUT ACTUAL TRIGGERS (ONLINE)
    % 1) Note down actual triggers - trig note
    AllSongsAllTrigsCollected(IFile).data_OnlineTrigs.ttimes=rd.ttimes;
    for iii=1:length(rd.pbname); % get notenum of triggers
        ind=findstr(rd.pbname{iii},'Templ');
        
        if isempty(ind) % then this is likley a manual trig. ignore it.
            assert(any(findstr(rd.pbname{iii},'Manual Trig')), 'PROBLEM, why si this not a real template or a manual trig?');

                    AllSongsAllTrigsCollected(IFile).data_OnlineTrigs.trignotes(iii)=nan;

                    continue
        end
        
        % == trig note
        AllSongsAllTrigsCollected(IFile).data_OnlineTrigs.trignotes(iii)=str2num(rd.pbname{iii}(ind+8));
        
    end
    
    
    % 2) Note down whether trigs were catch (either catch song or trial)
    AllSongsAllTrigsCollected(IFile).data_OnlineTrigs.catch_EitherSongOrTrial=rd.catch;
    AllSongsAllTrigsCollected(IFile).data_OnlineTrigs.IsCatchSong=rd.iscatch;
    % ================
    
    
    
    % 3) Finally, SAVE ACTUAL TRIGS (ONLINE) TO DIFFERENT FIELDNAMES IN RD FILE
    % Write Actual trigs (config FF) to rec file
    % NOTE: Used to use SimTrigInfoActual (and save this assuming this
    % is actual online trigs. BUT it is not necessarily accurate for 2
    % reasons: 1) uses config file FF thresh and 2) ignores Note Group
    % functionality - see above.  Therefore now saving online rd.ttimes
    % as the actual triggers.  It is accurate
    
    rd.ttimesActualTrig=rd.ttimes;
    
    for iii=1:length(rd.pbname); % get notenum of triggers
        ind=findstr(rd.pbname{iii},'Templ');
        
                if isempty(ind) % then this is likley a manual trig. ignore it.
            assert(any(findstr(rd.pbname{iii},'Manual Trig')), 'PROBLEM, why si this not a real template or a manual trig?');

        rd.trignoteActualTrig(iii)=nan;

                    continue
        end

        rd.trignoteActualTrig(iii)=str2num(rd.pbname{iii}(ind+8));
    end

%         rd.FreqValsActualTrig=SimTrigInfoActual(:,3); % left out.  can
%         get FF from comparing ttimes to offline ttimes


    
    % ================== Write All Offline trigs to rec file (i.e. all Detects)
    if ~isempty(SimTrigInfoAll); % i.e. this song contains at least one trig.
        [sortv,sorti]=sort(SimTrigInfoAll(:,1));
        SimTrigInfoAll=SimTrigInfoAll(sorti,:);
        rd.ttimes=SimTrigInfoAll(:,1);
        rd.trignote=SimTrigInfoAll(:,2);
        rd.FreqVals=SimTrigInfoAll(:,3);
        
    else
        rd.ttimes=[];
        rd.trignote=[];
        rd.FreqVals=[];
        
    end
    
        
    
    wrtrecf_LT_evtafv4(fn,rd,1); 
%     wrtrecf(fn,rd,1);
    

    % ========== old debug code
    %for iND=1:length(ND)
    %disp(['IND  = ',num2str(iND)]);
    %for ijk=1:length(ND(iND).TriggerTimes)
    %    disp(num2str(ND(iND).TriggerTimes(ijk),6));
    %end
    %disp([' ']);
    
    %rd=readrecf(fn);
    %tt=[rd.ttimes,zeros(size(rd.ttimes))];
    %for ijk=1:size(tt,1)
    %    pppp=findstr(rd.pbname{ijk},'Templ = ');
    %    strtemp=rd.pbname{ijk};
    %    tt(ijk,2)=str2num(strtemp(pppp+8:end));
    %end
    
    
    % Some Debug output i used early on
    %compare the sim trigs to the actual trigs
    %if (isempty(tt))
    %    disp(['No trigs in recfile : Nact = 0  - Nsim = ',num2str(length(ND(iND).TriggerTimes))]);
    %    continue;
    %end
    
    
    %pppp=find(tt(:,2)==iND-1);
    %if (length(pppp)==0)
    %    continue;
    %end
    %tttemp=tt(pppp,1);
    %if (length(tttemp)~=length(ND(iND).TriggerTimes))
    %    disp(['XTrig Times not the same size : Nact = ',num2str(length(tttemp)),' Nsim = ',num2str(length(ND(iND).TriggerTimes))]);
    %end
    %for ijk=1:length(ND(iND).TriggerTimes)
    %    [yy,iy]=min(abs(ND(iND).TriggerTimes(ijk)*1e3-tttemp));
    %    disp(['iND = ',num2str(iND)]);
    %    disp(['Act Trig Time = ',num2str(tttemp(iy)),' Sim Trigger Time = ',num2str(ND(iND).TriggerTimes(ijk)*1e3)]);
    %    disp(['Diff in ms = ',num2str( tttemp(iy) - ND(iND).TriggerTimes(ijk)*1e3)]);
    %end
end
end % loop over files




