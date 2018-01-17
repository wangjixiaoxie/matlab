function lt_batchsong_PlotSongTrigs(batchf, trignote, plotspec, plotAdjacSongIfTrig, ...
    plotPreviousSongIfTrig)
%% lt 1/12/18 (was script before) - plots song + trigger data from evtafv4 (+ online rec

%%

% plotspec = 1; % if 0, then plots sound dat, if 1, then spectrogram
% batchf = 'batchtmp';
% trignote = 0; % the DIO chan for the trigger (0 by default is WN)
% plotAdjacSongIfTrig=1; % if 1, then plots the song file directly after (searches directory), if there is trigger
% this is useful if expect trigger to last to next song file (e.g. magic
% glass experiment)

%% get list of all songs if needed
if plotAdjacSongIfTrig==1 | plotPreviousSongIfTrig==1
    % --- get list of songs in this dir
    
    Songsindir = dir('*.cbin');
    
    % -- only keep files (throw out dirs)
    Songsindir = Songsindir([Songsindir.isdir]==0);
    
    % -- sort by date
    [~, indtmp] = sort([Songsindir.datenum]);
    Songsindir = Songsindir(indtmp);
end


%% RUN
figcount=1;
subplotrows=4;
subplotcols=2;
fignums_alreadyused=[];
hfigs=[];

hsplots = [];

% ============================= go thur all songs in batch
fid = fopen(batchf);
songf = fgetl(fid);

while ischar(songf)
    
    % ----- start a new subplot
    [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
    title(songf);
    hsplots = [hsplots hsplot];
    
    % 1) ############################################# song
    [dat, Fs, DOFILT, ext]=ReadDataFile(songf,'0');
    x = [1:length(dat)]./Fs;
    
    % ----- plot dat or spectrogram
    if plotspec ==1
        lt_plot_spectrogram(dat, Fs, 0, 0);
    else
        plot(x, dat, '-k')
    end
    
    
    % 2) ######################################## trigger
    [dat, Fs, DOFILT, ext]=ReadDataFile(songf,'1');
    plot(x, dat, '-r')
    
    % -- ad trig time with rec file ...
    rd = readrecf_LT_evtafv4(songf);
    istrigMain =0;
    if ~isempty(rd.ttimes)
        for i=1:length(rd.ttimes)
            if rd.trignote(i)==trignote
                if rd.catch(i) ==1
                    plot(rd.ttimes(i)/1000, 0, '^g');
                else
                    plot(rd.ttimes(i)/1000, 0, '^r');
                    lt_plot_text(rd.ttimes(i)/1000 +0.05, 0, 'hit', 'r');
                    istrigMain =1;
                end
                
            end
            
        end
    end
    
    
    % 3) #####################################3 PLOT ADJACENT FILE?
    if istrigMain==1
        if plotAdjacSongIfTrig==1
            
            % --- locate this song in list of alls ongs
            indsong = find(strcmp({Songsindir.name}, songf));
            assert(length(indsong)==1, 'asdf');
            
            songnext = Songsindir(indsong+1).name;
            
            % ========================================= PLOT
            
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title(['[NEXTSONG] ' songnext]);
            hsplots = [hsplots hsplot];
            
            % 1) song
            [dat, Fs, DOFILT, ext]=ReadDataFile(songnext,'0');
            x = [1:length(dat)]./Fs;
            
            % ----- plot dat or spectrogram
            if plotspec ==1
                lt_plot_spectrogram(dat, Fs, 0, 0);
            else
                plot(x, dat, '-k')
            end
            
            % 2) trigger
            [dat, Fs, DOFILT, ext]=ReadDataFile(songnext,'1');
            plot(x, dat, '-r')
            
            % -- ad trig time with rec file ...
            rd = readrecf_LT_evtafv4(songnext);
            if ~isempty(rd.ttimes)
                for i=1:length(rd.ttimes)
                    if rd.trignote(i)==trignote
                        if rd.catch(i) ==1
                            plot(rd.ttimes(i)/1000, 0, '^g');
                        else
                            plot(rd.ttimes(i)/1000, 0, '^r');
                            lt_plot_text(rd.ttimes(i)/1000 +0.05, 0, 'hit', 'r');
                        end
                        
                    end
                    
                end
            end
        end
        
        if plotPreviousSongIfTrig==1
            
            % --- locate this song in list of alls ongs
            indsong = find(strcmp({Songsindir.name}, songf));
            assert(length(indsong)==1, 'asdf');
            
            songpre = Songsindir(indsong-1).name;
            
            % ========================================= PLOT
            
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title(['[PREVIOUSSONG] ' songpre]);
            hsplots = [hsplots hsplot];
            
            % 1) song
            [dat, Fs, DOFILT, ext]=ReadDataFile(songpre,'0');
            x = [1:length(dat)]./Fs;
            
            % ----- plot dat or spectrogram
            if plotspec ==1
                lt_plot_spectrogram(dat, Fs, 0, 0);
            else
                plot(x, dat, '-k')
            end
            
            % 2) trigger
            [dat, Fs, DOFILT, ext]=ReadDataFile(songpre,'1');
            plot(x, dat, '-r')
            
            % -- ad trig time with rec file ...
            rd = readrecf_LT_evtafv4(songpre);
            if ~isempty(rd.ttimes)
                for i=1:length(rd.ttimes)
                    if rd.trignote(i)==trignote
                        if rd.catch(i) ==1
                            plot(rd.ttimes(i)/1000, 0, '^g');
                        else
                            plot(rd.ttimes(i)/1000, 0, '^r');
                            lt_plot_text(rd.ttimes(i)/1000 +0.05, 0, 'hit', 'r');
                        end
                        
                    end
                    
                end
            end
            
        end
    end
    
    % ------------ next song
    songf = fgetl(fid);
end

linkaxes(hsplots, 'y');
