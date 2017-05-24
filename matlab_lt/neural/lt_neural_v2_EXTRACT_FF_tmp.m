function lt_neural_v2_EXTRACT_FF(SummaryStruct, FFparamsAll, overWrite, ...
    plotSpec, plotOnSong, plotSyl)

%%

% overWrite = 1;
% plotSpec = 0; % to plot raw spec overlayed with PC and windows.
% plotOnSong = 2; % will only start plotting spec once hit this song num.
% plotSyl = ''; % to focus on just one syl. NOT DONE YET

if overWrite == 1
    assert(strcmp(input('sure you want to overwrite? (y or n) ', 's'), 'y'), 'stopping!!');
end

%% lt 3/22/17 - extract FF for all neurons in Summary struct. make note in struct.
prepad=0.015; postpad=0.015; % get 15ms pre and post (acoustic dat) [DO NOT CHANGE!! - since
% timebin windows are defined relative to this onset

filenamedatFF = 'extractFF';
filenamedatPC = 'extractPC';
filenameparams = 'extractFF_params';

%%

Numbirds = length(SummaryStruct.birds);
for z = 1:Numbirds
    Numneurons = length(SummaryStruct.birds(z).neurons);
    birdname = SummaryStruct.birds(z).birdname;
    
    for zz = 1:Numneurons
        
        datstruct = SummaryStruct.birds(z).neurons(zz);
        cd(datstruct.dirname);
        
        % ==================== EXTRACT
        % load dat for neuron
        cd ..
        batchf = datstruct.batchfilename;
        chan = datstruct.channel;
        [SongDat, NeurDat, ~] = lt_neural_ExtractDat(batchf, chan, 0);
        
        % 1) First check whether FF has already been extracted
        previouslydone = exist([filenamedatFF '.mat'] , 'file'); % 2 if done, 0 if not
        
        
        if previouslydone==2 % decide whether to do
            
            % ============================== HAS PARAMS CHANGED?
            paramsSame=1;
            tmp = load('extractFF_params');
            
            if  ~strcmp(SongDat.AllLabels, tmp.Params.allLabels)==1 | ...
                    ~all(SongDat.AllOnsets == tmp.Params.AllOnsets)
                paramsSame=0;
            end
            
            % === ff params
            ind = find(strcmp({FFparamsAll.bird.birdname}, birdname)); % find info for this syl
            cell_of_freqwinds=FFparamsAll.bird(ind).FFparams.cell_of_freqwinds; % then is not learning
            
            % if is learning experiement, then use diff, since WN will be over
            % a syl
            if isempty(datstruct.LEARN_WNonDatestr)
                cell_of_FFtimebins=FFparamsAll.bird(ind).FFparams.cell_of_FFtimebins;
            else
                cell_of_FFtimebins=FFparamsAll.bird(ind).FFparams.cell_of_FFtimebins_DurLearn; % then is not learning
            end
            
            if ~isequal(tmp.Params.cell_of_freqwinds, cell_of_freqwinds) | ...
                    ~isequal(tmp.Params.cell_of_FFtimebins, cell_of_FFtimebins)
                paramsSame=0;
            end
            % ==================================================================
            
            
            % ----
            if overWrite == 0 & paramsSame ==1
                disp(['Skipping bird ' num2str(z) ' neuron ' num2str(zz) '[PREVIOUS DAT FOUND]']);
                % skip
                continue
            end
            
            if overWrite == 0 & paramsSame ==0
                disp(['Overwriting bird ' num2str(z) ' neuron ' num2str(zz) '[since params changed]']);
            end
            
            % 2) If not extracted or told to overwrite, reextract
            if overWrite ==1
                disp(['Overwriting bird ' num2str(z) ' neuron ' num2str(zz) '[PREVIOUS DAT FOUND]']);
            end
        else
            disp(['Extracting for the first time:  bird ' num2str(z) ' neuron ' num2str(zz)]);
        end
        
        
        %% ============= RUN
        % load dat for neuron
        cd ..
        batchf = datstruct.batchfilename;
        chan = datstruct.channel;
        [SongDat, NeurDat, ~] = lt_neural_ExtractDat(batchf, chan, 1);

        % -- extract infor for this bird.
        ind = find(strcmp({FFparamsAll.bird.birdname}, birdname)); % find info for this syl
        cell_of_freqwinds=FFparamsAll.bird(ind).FFparams.cell_of_freqwinds; % then is not learning
        
        % if is learning experiement, then use diff, since WN will be over
        % a syl
        if isempty(datstruct.LEARN_WNonDatestr)
            cell_of_FFtimebins=FFparamsAll.bird(ind).FFparams.cell_of_FFtimebins;
        else
            cell_of_FFtimebins=FFparamsAll.bird(ind).FFparams.cell_of_FFtimebins_DurLearn; % then is not learning
        end
        
        % ========= Go thru all labeled syls, for each, collect data and store
        % in output data structure
        numsyls=length(SongDat.AllLabels);
        FFvals = NaN(1, numsyls, 'single');
        PCvals = cell(1, numsyls);
        Tbasevals = cell(1, numsyls);
        fs=NeurDat.metaDat(1).fs;
        
        for j=1:numsyls
            syl=SongDat.AllLabels(j);
            
            if strcmp(syl, '-') % only collect labeled syls
                continue
            end
            
            % --- acoustic data
            onset_sec=SongDat.AllOnsets(j);
            onset_samp=int64(onset_sec*fs);
            
            offset_sec=SongDat.AllOffsets(j);
            offset_samp=int64(offset_sec*fs);
            
            syldat=SongDat.AllSongs(onset_samp-prepad*fs:offset_samp+postpad*fs);
            
            % Pitch contour and FF ----------------------
            ind=find(strcmp(syl, cell_of_freqwinds));
            collectFF=1;
            if isempty(ind)
                % tel;l use
                %             disp(['NO FREQ WINDOW SPECIFIED FOR SYL ' syl ' !!']);
                collectFF=0;
            elseif isempty(cell_of_freqwinds{ind+1})
                collectFF=0;
            end
            
            % -- defaults, sets output (e.g., ff) to this if don't collect
            % actual FF
            if collectFF==1
                ind=find(strcmp(syl, cell_of_freqwinds));
                F_high=cell_of_freqwinds{ind+1}(2);
                F_low=cell_of_freqwinds{ind+1}(1);
                
                ind=find(strcmp(cell_of_FFtimebins, syl));
                mintime=cell_of_FFtimebins{ind+1}(1); % sec
                maxtime=cell_of_FFtimebins{ind+1}(2);
                
                currsong = SongDat.AllSongNum(j);
                if plotSpec==1
                    % only plot if reached desired song num
                    if currsong<plotOnSong
                        [FF, PC, T]= lt_calc_FF(syldat, fs, [F_low F_high], [mintime maxtime], 0);
                    else
                        disp(syl);
                        [FF, PC, T]= lt_calc_FF(syldat, fs, [F_low F_high], [mintime maxtime], 1);
                    end
                else
                    [FF, PC, T]= lt_calc_FF(syldat, fs, [F_low F_high], [mintime maxtime], 0);
                end
                
                
                % ---------- OUTPUT
                FFvals(j) = FF;
                PCvals{j} = PC;
                Tbasevals{j} = T;
                %             plot(T, PC, 'r', 'LineWidth', 2);
            end
            
            
            
            % ---- What song file this came from?
            %             globalOnsetTime=SongDat.AllOnsets(j); % sec
            %             globalOnsetSamp=globalOnsetTime*fs;
            %             cumOnsSamps=cumsum([0 NeurDat.metaDat.numSamps]);
            %             songind=find((globalOnsetSamp-cumOnsSamps)>0, 1, 'last');
            %             songfname=NeurDat.metaDat(songind).filename;
        end
        
        % ------- SAVE
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


