function lt_disp_random_songs_syls(batch,SongOrSyl, israndom, RandNum, plotFB, FBnote, FBdur, varargin)
%% LT 12/20/14 - Displays random (or all) songs, or random (or all) of a designated syllable
% displaying spectrograms same way as evsonganaly. modified from that code.

% Example inputs:
% batch='batch.labeled.all';
% SongOrSyl=1; if 1, then plots syllables (then needs varargin, see below).
%   If 0, then plots entire songs. Don't need varargin
% israndom=1; % if 1, then takes random number of songs/syls to plot
% RandNum=10; % if random, then specificy how to look at?
% plotFB=1; % if 1, then plots triggers (actual, from rec file), then
% FBnote=0; % specificy the note of FB (i.e. evtafv4 note);
% FBdur=0.2; % in sec, duration of FB, for plotting purposes. Use 0 to plot
%   just a dot

% varargin is if plotting syllable, see above
% SylTarg='b'; varargin{1}
% PreNote='bcc'; varargin{2}
% PreDur=0.6; % in sec, to capture "motif", change durations; varargin{3}
% PostSylOnsetDur=0.2; varargin{4}

% EXAMPLES:
% lt_disp_random_songs_syls(batch,0,1,10); % to show 10 random songs
%
% lt_disp_random_songs_syls(batch,1,1,30,1,0,0,'b','dcc',0.2,0.2); % to
% show random syls
%
% lt_disp_random_songs_syls(batch,1,0,0,1,0,0,'b','dcc',0.2,0.2); % to show
% all syls


if size(varargin,2)>0;
    SylTarg=varargin{1};
    PreNote=varargin{2};
    PreDur=varargin{3};
    PostSylOnsetDur=varargin{4}; % in sec, to capture "motif", change durations
end

Fs=32000;


    % FILL IN NONEXISTENT ARGS
    if ~exist('plotFB','var');
        plotFB=0;
%         FBnote=0;
%         FBdur=0;
    end
    
    if ~exist('PreDur','var');
        PreDur=0.016;
    end
    
      if ~exist('PostSylOnsetDur','var');
        PostSylOnsetDur=0.02;
    end
      
        
if SongOrSyl==0; % then plot songs
    
    
    
    
    
    
    %% EXTRACT FILENAMES FROM BATCH and put into cell array
    
    fid1=fopen(batch);
    fn=fgetl(fid1);
    i=1;
    
    while ischar(fn);
        fnames{i}=fn; % put filenames into cell array
        
        i=i+1;
        fn=fgetl(fid1);
    end
    
    
    % GET RANDOM SUBSET IF DESIRED
    if israndom==1;
        if length(fnames)>RandNum; % only get random subset if enough songs
            TT=randperm(length(fnames)); % perm list of indices
            RandInds=TT(1:RandNum); % pick certain number of ran inds.
            RandInds=sort(RandInds); % get in chrono order
            fnames=fnames(RandInds); % extract only those random fnames
        end
    end
    
    
    
    %% LOAD AND PLOT ALL SONGS
    NumSongs=length(fnames);
    
    % plot parameters
    MaxPlotsPerFig=8;
    row_plots=2;
    col_plots=4;
    num_figures=ceil(NumSongs/MaxPlotsPerFig);
    
    % WARNING - if too many figures
    if num_figures>10;
        if strcmp(input(['Will create ' num2str(num_figures) ' figures. Type y to continue'],'s'),'y')==0;
            sacasdfasd; % to cause error in code.
        end
    end
    
    % Initiate figures
    for i=1:num_figures;
        
        hfig1(i)=figure;
    end
    
    
    c=1;
    for i=1:length(fnames);
        
        % LOAD song file
        fname=fnames{i};
        [dat, Fs, DOFILT, ext]=ReadDataFile(fname,'0');
        
        
        % SMOOTH sound file
        [sm,sp,t,f]=SmoothData(dat,Fs,DOFILT,'‘hanningfirff’'); % ends up doing buttor, with filtfilt.
        
        % REMOVE unwanted frequencies and times
        sp=sp(10:140,:);
        f=f(10:140);
        
        sp=sp(:,312:end); % remove 1st sec.
        t=t(312:end);
        
        
        
        % TAKE LOG of sp
        % first, convert any sp values of 0 to non-zero(to the lowest value present);
        % solves problem of taking log of 0
        pp=find(sp>0);
        mntmp = min(min(sp(pp)));
        pp=find(sp==0);
        sp(pp) = mntmp;
        
        % second, take log
        sptemp=log(sp);
        sptemp = sptemp - min(min(sptemp));
        sptemp = uint8(((2^8) - 1)*(sptemp./max(max(sptemp)))); % SAVE SOME MEMORY 8X less than 64 bit double
        
        
        % PLOT
        % Prepare subplot slot
        fignum=ceil(c/MaxPlotsPerFig);
        subplotnum=mod(c,MaxPlotsPerFig);
        if subplotnum==0;
            subplotnum=MaxPlotsPerFig;
        end
        
        figure(hfig1(fignum));
        lt_plot_format;
        h1(c)=subplot(col_plots,row_plots,subplotnum); hold on;
        
        % Plot
        imagesc(t, f, sptemp);
        
        title(fname);
        axis([t(1) t(end) f(1) f(end)]);
        
        
        % MARK AS NOTCATCH
        try
        rd=readrecf(fname);
        if rd.iscatch==0;
            rectangle('Position',[t(1)+0.07, f(1)+70, t(end)-t(1)-0.14, f(end)-f(1)-140],'LineWidth',2, 'EdgeColor','r')
        end
        catch
            % i.e. rec file doesnt exist
        end
        
        % ANNOTATE WITH FB TIMES
        % note: is assuming evtaf v4, so has multiple notes.
        if plotFB==1
            for jj=1:length(rd.ttimes);
                TrigTimeSec=rd.ttimes(jj)/1000;
                TrigTimeSec=TrigTimeSec-0.012; % empircally chosen adjustment
                if FBdur==0; % i.e. plot a dot
                    plot(TrigTimeSec,f(1),'^b','MarkerFaceColor','b','MarkerSize',5);
                else
                    plot(TrigTimeSec,f(1),'^b','MarkerFaceColor','b','MarkerSize',5);
                    line([TrigTimeSec TrigTimeSec+FBdur], [f(1)+150 f(1)+150],'LineWidth',3);
                end
            end
        end
        
        
        % Advance by 1
        c=c+1;
        
        if fignum > 10;
           if strcmp(input('enter "c" to close current figs ', 's'), 'c')
               close all;
               c = 1;
           end
        end
        
    end
    
    linkaxes(h1,'y')
    
end

%% TO PLOT SYLLABLE
% can choose to plot by syllable or by stim presence

if SongOrSyl==1;
    
    PreTime=-PreDur; % the 0.016 is becuase the 1st window give a value at the middle of the window. The window is 32ms, so it's value is 16ms later than desired.
    frequency_range=''; % not calculating FF
    DatDur=(PreDur+PostSylOnsetDur)*Fs; % how much data (time wise) to get?
    
    
    % GET RAW SYL DATA
    SylData=findwnote2tw_v4_LT(batch, SylTarg,PreNote,PreTime,...
        frequency_range,DatDur,1,1,'obs0',0);
    
    SylNum=size(SylData,2);
    
    
    % GET RANDOM SUBSET IF DESIRED
    if israndom==1;
        if SylNum>RandNum;
            TT=randperm(SylNum);
            RandInds=TT(1:RandNum);
            RandInds=sort(RandInds); % get chrono order
            SylData=SylData(RandInds);
            
            SylNum=RandNum;
        end
        
        
    end
    
    
    % plot parameters
    MaxPlotsPerFig=25;
    row_plots=ceil(sqrt(MaxPlotsPerFig));
    col_plots=row_plots;
    num_figures=ceil(SylNum/MaxPlotsPerFig);
    
    
    % WARNING - too many figures?
    if num_figures>10;
        if strcmp(input(['Will create ' num2str(num_figures) ' figures. Type y to continue'],'s'),'y')==0;
            sacasdfasd; % to cause error in code.
        end
    end
    
    
    % Initiate figures
    for i=1:num_figures;
        
        hfig2(i)=figure;
    end
    
    
    % PLOT DATA
    
    c=1;
    for i=1:SylNum;
        dat=SylData(i).datt;
        % get filename without notmat
        xx=strfind(SylData(i).fn,'.cbin.not.mat');
        fname=[SylData(i).fn(1:xx-1) '; Note ' num2str(SylData(i).NotePos)];
        
        
        % SMOOTH sound file
        [sm,sp,t,f]=SmoothData(dat,Fs,1,'‘hanningfirff’'); % ends up doing buttor, with filtfilt.
        
        % REMOVE unwanted frequencies and times
        sp=sp(10:140,:);
        f=f(10:140);
        
        % sp=sp(:,312:end); % remove 1st sec.
        % t=t(312:end);
        
        
        
        % TAKE LOG of sp
        % first, convert any sp values of 0 to non-zero(to the lowest value present);
        % solves problem of taking log of 0
        pp=find(sp>0);
        mntmp = min(min(sp(pp)));
        pp=find(sp==0);
        sp(pp) = mntmp;
        
        % second, take log
        sptemp=log(sp);
        sptemp = sptemp - min(min(sptemp));
        sptemp = uint8(((2^8) - 1)*(sptemp./max(max(sptemp)))); % SAVE SOME MEMORY 8X less than 64 bit double
        
        
        % PLOT
        % Prepare subplot slot
        fignum=ceil(c/MaxPlotsPerFig);
subplotnum=mod(c,MaxPlotsPerFig);
        if subplotnum==0;
            subplotnum=MaxPlotsPerFig;
        end
        
        figure(hfig2(fignum));
        lt_plot_format;
        h2(c)=subplot(col_plots,row_plots,subplotnum); hold on;
        
        % Plot
        imagesc(t, f, sptemp);
        
        title(fname);
        axis([t(1) t(end) f(1) f(end)]);
        
        
        % MARK AS NOTCATCH
        % if rd.iscatch==0;
        %     rectangle('Position',[t(1)+0.07, f(1)+70, t(end)-t(1)-0.14, f(end)-f(1)-140],'LineWidth',2, 'EdgeColor','r')
        % end
        
        
        % ANNOTATE WITH NOTE TIMES
        % note: is assuming evtaf v4, so has multiple notes.
        plot(PreDur,f(1),'ok','MarkerFaceColor','k','MarkerSize',5);
        
        
        % PLOT TRIGGER IF ANY OCCURED
        if plotFB==1;
            rd=readrecf(SylData(i).fn);
            
            % what is window of data plotted?
            SylOns_sec=SylData(i).ons(SylData(i).NotePos)/1000; % onset of target syl, in sec
            
            t1=SylOns_sec-PreDur; % onset of window
            t2=SylOns_sec+PostSylOnsetDur; % offset
            
            % find triggers within that window
            TrigTimes_sec=rd.ttimes/1000;
            TrigsToPlot=find(TrigTimes_sec>t1 & TrigTimes_sec<t2);
            
            if ~isempty(TrigsToPlot)>0;
                for jj=1:length(TrigsToPlot);
                    TrigOnsRelEpoch=rd.ttimes(TrigsToPlot(jj))/1000-t1;
                    TrigOnsRelEpoch=TrigOnsRelEpoch-0.012; % empirically determined adjustment
                    if FBdur==0; % i.e. plot a dot
                        plot(TrigOnsRelEpoch,f(1),'^r','MarkerFaceColor','r','MarkerSize',5);
                    else
                        plot(TrigOnsRelEpoch,f(1),'^r','MarkerFaceColor','r','MarkerSize',5);
                        line([TrigOnsRelEpoch TrigOnsRelEpoch+FBdur], [f(1)+150 f(1)+150],'LineWidth',3, 'Color','r');
                    end
                end
            end
        end
        
        
        % Advance by 1
        c=c+1;
    end
    linkaxes(h2,'xy')
    
    
end


%% TO PLOT STIMS



















