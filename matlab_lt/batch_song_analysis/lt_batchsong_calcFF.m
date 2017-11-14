function lt_batchsong_calcFF(ListOfDirs, ListOfBatch, FFparams, plotAllPC, plotEachSyl);

%% lt 11/13/17 - given batch and labeled data, calculates FF for all syls.
% saves in a mat file in same dir as song.
% works for both cbin and intan files

% NOTE: only labels lowercase syllables.

%% === params

% ListOfDirs = {...
%     '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/111317_Afternoon_DirUndir/UNDIR', ...
%     '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/111317_Afternoon_DirUndir/DIR'};
% ListOfBatch = {...
%     'batchall', ...
%     'batchall'};
%
% FFparams.cell_of_freqwinds={'h', [2340 3160], 'b', [2340 3160]}; % 'j', [950 1450], 'l', [1200 1600], 't', [3590 4960]
% FFparams.cell_of_FFtimebins={'h', [0.043 0.058], 'b', [0.053 0.07]}; % 'j', [0.04 0.045], 'l', [0.035 0.039], 't', [0.026 0.033], ...
% FFparams.cell_of_FFtimebins_DurLearn={'h', [0.034 0.038], 'b', [0.053 0.07]}; % WN on H
%
% plotAllPC =1;
% plotEachSyl = 1;

%% defaults

paddur = 0.015; % in sec DONT CHANGE!!

assert(paddur==0.015, 'oproblem, timing messed up');


%% ==== desired syls
sylstoget = unique([FFparams.cell_of_FFtimebins{1:2:end}]);

%% run

if plotAllPC ==1
    
    % will collect all PC
    PCstruct = struct;
    
    for i=1:length(sylstoget)
        PCstruct.(sylstoget(i)).PC = {};
        PCstruct.(sylstoget(i)).T = {};
        PCstruct.(sylstoget(i)).ff = [];
        
    end
end


for i=1:length(ListOfDirs)
    
    dirname = ListOfDirs{i};
    batchf = ListOfBatch{i};
    
    cd(dirname);
    
    % --- go thru all songs in batchf
    fid = fopen(batchf);
    fname = fgetl(fid);
    
    while ischar(fname)
        disp(fname);
        
        % ------ if previously done, overwrite?
        % TO DO!!!!!!!!!!!!!!
        
        % --- if no notmat, then skip
        if ~exist([fname '.not.mat'], 'file')
            fname = fgetl(fid);
            disp('==== SKIP - no notmat')
            continue
        end
        
        tmp = load([fname '.not.mat']);
        inds = regexp(tmp.labels, ['[' sylstoget ']']);
        
        if isempty(inds)
            % then desired syl is not labeled, skip
            fname = fgetl(fid);
            disp(' ====== SKIP - desired syl not labeled');
            continue
        end
        
        % ------- LOAD SONG FILE
        [frequency_parameters, songdat] = pj_readIntanNoGui_AudioOnly(fname);
        fs = frequency_parameters.amplifier_sample_rate;
        
        % =============== EXTRACT FF
        FFall = nan(1,length(tmp.labels));
        
        % --- for plotting each syl
        figcount=1;
        subplotrows=6;
        subplotcols=4;
        fignums_alreadyused=[];
        hfigs=[];
        
        
        
        for j=inds
            
            % =============== extract FF
            syl = tmp.labels(j);
            disp(syl);
            
            % ---- EXTRACT SOUND SEGEMENT
            ons = tmp.onsets(j) - paddur*1000;
            off = tmp.offsets(j) + paddur*1000;
            
            ons = round(fs*ons/1000); % - convert to samps
            off = round(fs*off/1000);
            
            syldat = songdat(ons:off);
            
            
            % ---- EXTRACT PARAMS
            ind=find(strcmp(syl, FFparams.cell_of_freqwinds));
            F_high=FFparams.cell_of_freqwinds{ind+1}(2);
            F_low=FFparams.cell_of_freqwinds{ind+1}(1);
            
            ind=find(strcmp(FFparams.cell_of_FFtimebins, syl));
            mintime=FFparams.cell_of_FFtimebins{ind+1}(1); % sec
            maxtime=FFparams.cell_of_FFtimebins{ind+1}(2);
            
            
            % ---- RUN, EXTRACT FF
            [ff, PC, T]= lt_calc_FF(syldat, fs, [F_low F_high], [mintime maxtime], 0);
            
            % =================== OUTPUT
            FFall(j) = ff;
            
            if plotAllPC==1
                PCstruct.(syl).PC = [PCstruct.(syl).PC PC];
                PCstruct.(syl).T = [PCstruct.(syl).T T];
                PCstruct.(syl).ff = [PCstruct.(syl).ff ff];
            end
            
            if plotEachSyl==1
                [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
                title([syl ', pos:' num2str(j)]);
                lt_plot_spectrogram(syldat, fs, 0, 0, [], [], 0);
                plot(T, PC, 'w', 'LineWidth', 2);
                
                % --- overlay FF and time windows
                line([mintime mintime], ylim, 'Color', 'k');
                line([maxtime maxtime], ylim, 'Color', 'k');
                
                line(xlim, [F_high F_high], 'Color', 'k');
                line(xlim, [F_low F_low], 'Color','k');
                
                line([mintime maxtime], [ff ff], 'Color', 'k', 'LineWidth', 2);
            end
        end
        
        if plotEachSyl==1
            lt_subtitle(fname);
            pause;
            close all;
        end
        
        
        % ======================= SAVE FF IN A MAT FILE FOR THIS SONG
        tstamp = lt_get_timestamp(0);
        
        FFstruct = struct;
        FFstruct.FFall = FFall;
        FFstruct.FFparams = FFparams;
        FFstruct.extraction_tstamp = tstamp;
        
        fnamesave = [fname '.calcff.mat'];
        save(fnamesave, 'FFstruct');
        
        fname = fgetl(fid);
    end
end

%% ====== plot all PCs if wanted
if plotAllPC==1
    
    figcount=1;
    subplotrows=6;
    subplotcols=3;
    fignums_alreadyused=[];
    hfigs=[];
    
    
    for i=1:length(sylstoget)
        
        syl = sylstoget(i);
        
        % ---- EXTRACT PARAMS
        ind=find(strcmp(syl, FFparams.cell_of_freqwinds));
        F_high=FFparams.cell_of_freqwinds{ind+1}(2);
        F_low=FFparams.cell_of_freqwinds{ind+1}(1);
        
        ind=find(strcmp(FFparams.cell_of_FFtimebins, syl));
        mintime=FFparams.cell_of_FFtimebins{ind+1}(1); % sec
        maxtime=FFparams.cell_of_FFtimebins{ind+1}(2);
        
        % ---
        PCcell = PCstruct.(syl).PC;
        Tcell = PCstruct.(syl).T;
        [fignums_alreadyused, hfigs, figcount, hsplot]=lt_plot_MultSubplotsFigs('', subplotrows, subplotcols, fignums_alreadyused, hfigs, figcount);
        title(syl);
        
        for j=1:length(PCcell)
            pc = PCcell{j};
            t = Tcell{j};
            ff = PCstruct.(syl).ff(j);
            
            plot(t, pc);
            
            % --- overlay FF and time windows
            line([mintime mintime], ylim, 'Color', 'r');
            line([maxtime maxtime], ylim, 'Color', 'r');
            
            line(xlim, [F_high F_high], 'Color', 'r');
            line(xlim, [F_low F_low], 'Color','r');
            
            line([mintime maxtime], [ff ff], 'Color', 'm');
            
            
        end
    end
end
