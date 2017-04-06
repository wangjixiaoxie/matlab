function lt_neural_v2_EXTRACT_WNhit(SummaryStruct, FFparamsAll, overWrite, ...
    plotSpec, plotOnSong, plotSyl)
%% LT modified from extractFF - here extracts whether hit by WN for each syl
% INCOMPLETE!!!! - stopped at line 64, was trying to decide how muhc of syl
% data to take to search for WN (e.g. just syl, or with some pre and post).
% Decided will just do this during analysis - as then know time of Wn
% relative to amount of extracted data.

refractPeriod = 0.05; % any diffs between max audio datapoints considered same unless diff greater than this.
%% 

% overWrite = 1;
% plotSpec = 0; % to plot raw spec overlayed with PC and windows.
% plotOnSong = 2; % will only start plotting spec once hit this song num.
% plotSyl = ''; % to focus on just one syl. NOT DONE YET

%% lt 3/22/17 - extract FF for all neurons in Summary struct. make note in struct.
prepad=0.015; postpad=0.015; % get 15ms pre and post (acoustic dat) [DO NOT CHANGE!! - since
% timebin windows are defined relative to this onset

filenamedatFF = 'extractWNhit';
filenameparams = 'extractWNhit_params';

%%

Numbirds = length(SummaryStruct.birds);
for z = 1:Numbirds
    Numneurons = length(SummaryStruct.birds(z).neurons);
    birdname = SummaryStruct.birds(z).birdname;
    
    for zz = 1:Numneurons
        
        datstruct = SummaryStruct.birds(z).neurons(zz);
        
        % ==================== EXTRACT
        cd(datstruct.dirname);
        
        % 1) First check whether FF has already been extracted
        previouslydone = exist([filenamedatFF '.mat'] , 'file'); % 2 if done, 0 if not
        
        if overWrite == 0 & previouslydone == 2
            disp(['Skipping bird ' num2str(z) ' neuron ' num2str(zz) '[PREVIOUS DAT FOUND]']);
            % skip
            continue
        end
        
        % 2) If not extracted or told to overwrite, reextract
        if overWrite ==1 & previouslydone ==2
            disp(['Overwriting bird ' num2str(z) ' neuron ' num2str(zz) '[PREVIOUS DAT FOUND]']);
        end
        if previouslydone ==0
            disp(['Extracting for the first time:  bird ' num2str(z) ' neuron ' num2str(zz)]);
        end
        
        
        % load dat for neuron
        cd ..
        batchf = datstruct.batchfilename;
        chan = datstruct.channel;
        extractSound=1;
        [SongDat, NeurDat, ~] = lt_neural_ExtractDat(batchf, chan, extractSound);
          
% === TAKE ENTIRE SONG CONCATED FILE, AND EXTRACT TIMES OF ALL WN
        tmp = find(SongDat.AllSongs > 3.2);


        %% ------- SAVE
        cd(datstruct.dirname);
        
        % - dat
        save(filenamedatFF, 'FFvals');
        save(filenamedatPC, 'PCvals');
 
      % - params
        Params = struct;
        Params.cell_of_FFtimebins = cell_of_FFtimebins;
        Params.cell_of_freqwinds = cell_of_freqwinds;
        Params.allLabels = SongDat.AllLabels;
        Params.AllOnsets = SongDat.AllOnsets;
        Params.AllOffsets = SongDat.AllOffsets;
        
        Params.tbin = T(2)-T(1);
        Params.tfirstbin = T(1);
        save(filenameparams, 'Params');
        
    end
end  

%% troubleshooting
if (0)
    %% load and plot all PCs for all syls for this neuron
    % go to save folder and run (once for each neuron)
    figcount=1;
    subplotrows=5;
    subplotcols=3;
    fignums_alreadyused=[];
    hfigs=[];
    
    load extractPC;
    load extractFF_params
    
    syltypes = unique(Params.allLabels);
    numsyls =length(unique(Params.allLabels));
    
    for j=1:numsyls;
        
        sylname = syltypes(j);
        
        
        ind = strfind(Params.allLabels, sylname);
        
        if ~isempty(ind)
            
            [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
            title(sylname);
            for i=1:length(ind)
                
                pc = PCvals{ind(i)};
                t = Params.tfirstbin:Params.tbin:Params.tfirstbin+Params.tbin*(length(pc)-1);
                plot(t, pc, ':');
                
            end
            try
            ind = find(strcmp(Params.cell_of_freqwinds, sylname));
            freqmin = Params.cell_of_freqwinds{ind+1}(1);
            freqmax = Params.cell_of_freqwinds{ind+1}(2);
            line(xlim, [freqmin freqmin]);
            line(xlim, [freqmax freqmax]);
            
            ind = find(strcmp(Params.cell_of_FFtimebins, sylname));
            tmin = Params.cell_of_FFtimebins{ind+1}(1);
            tmax = Params.cell_of_FFtimebins{ind+1}(2);
            line([tmin tmin], ylim);
            line([tmax tmax], ylim);
            catch err
            end
        end
    end
end
    

