%% calc FF for all syls, save next to song file
clear all; close all;

ListOfDirs = {...
    '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/111317_Afternoon_DirUndir/UNDIR', ...
    '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/111317_Afternoon_DirUndir/DIR'};
ListOfBatch = {...
    'batchall', ...
    'batchall'};

FFparams.cell_of_freqwinds={'h', [3000 3800], 'b', [2900 3800]}; % 'j', [950 1450], 'l', [1200 1600], 't', [3590 4960]
FFparams.cell_of_FFtimebins={'h', [0.033 0.042], 'b', [0.033 0.045]}; % 'j', [0.04 0.045], 'l', [0.035 0.039], 't', [0.026 0.033], ...

plotAllPC =1;
plotEachSyl = 0;

lt_batchsong_calcFF(ListOfDirs, ListOfBatch, FFparams, plotAllPC, plotEachSyl);



%% = extract FF and Time vectors from .calcFF.mat
% ================== METHOD 1) EXTRACTS ALL SINGLE SYLS

ListOfDirs = {...
    '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/111317_Afternoon_DirUndir/UNDIR', ...
    '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/111317_Afternoon_DirUndir/DIR'};
ListOfBatch = {...
    'batchall', ...
    'batchall'};

%%
DATSTRUCT = struct;
count=1;
for i=1:length(ListOfBatch)
    
    dirname = ListOfDirs{i};
    batchf = ListOfBatch{i};
    
    cd(dirname);
    
    % --- go thru all songs in batchf
    fid = fopen(batchf);
    fname = fgetl(fid);
    
    while ischar(fname)
        disp(fname);
        
        % ================ check if calFF.mat exists
        if ~exist([fname '.calcff.mat'], 'file')
            fname = fgetl(fid);
            disp('==== SKIP - no notmat')
            continue
        end
        
        % ============= extract FF and time for all syls
        notmat = load([fname '.not.mat']);
        calcff = load([fname '.calcff.mat']);
        
        if any(strfind(fname, '.rhd'));
            [dtnum datestring]=lt_neural_fn2datenum(fname);
        else
            disp('PROBLEM!!! - name extraction, do for evtaf');
        end
        
        % ============ iterate over all syls
        for j=1:length(notmat.labels)
            
            % --- skip if is missing ff val
            if isnan(calcff.FFstruct.FFall(j))
                %                 disp(['skipped ' notmat.labels(j)]);
                continue
            end
            
            % ============= output
            DATSTRUCT.rendnum(count).ff = calcff.FFstruct.FFall(j);
            DATSTRUCT.rendnum(count).syl = notmat.labels(j);
            DATSTRUCT.rendnum(count).fname = fname;
            DATSTRUCT.rendnum(count).dirname = dirname;
            DATSTRUCT.rendnum(count).time_song_SecRes = dtnum;
            DATSTRUCT.rendnum(count).time_withinsong = notmat.onsets(j)/1000;
            
            % ------
            count = count+1;
        end
        
        fname = fgetl(fid);
    end
end

%% === METHOD 2 - extracts specific motifs of interest
clear all; close all;
ListOfDirs_UNDIR = { ...
    '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/111317_Afternoon_DirUndir/UNDIR'};
ListOfDirs_DIR = {...
    '/bluejay5/lucas/birds/pu69wh78/NEURAL/110917_RALMANOvernightLearn1/111317_Afternoon_DirUndir/DIR'};

ListOfBatch = {...
    'batchall', ...
    'batchall'}; % should be in order of UNDIR --> DIR

MotifsToExtract = {'ab(h)', 'jb(h)'};

DATSTRUCT = lt_batchsong_extractFF(ListOfDirs_UNDIR, ListOfDirs_DIR, ListOfBatch, MotifsToExtract);

%% ==================== PLOT
for i = 1:length(MotifsToExtract)
    
    ffvals = [DATSTRUCT.motif(i).rendnum.ff];
    tvals = [DATSTRUCT.motif(i).rendnum.datenum_song_SecRes];
    isDir = [DATSTRUCT.motif(i).rendnum.isDIR];
    
    % -- convert tvals to days from start
    firstday = datestr(floor(tvals), 'ddmmmyyyy');
    tvals = lt_convert_EventTimes_to_RelTimes(firstday, tvals);
    tvals = tvals.FinalValue;
    
    % ---- UNDIR
    figure; hold on; 
    plot(tvals(isDir==0), ffvals(isDir==0), 'ok');
    
    % --- DIR
    plot(tvals(isDir==1), ffvals(isDir==1), 'ob');

    title(MotifsToExtract{i});
end
    
    
    
    
    
    
    
    
    
    
    

