function lt_neural_concatOneChan(batchf, channel_board)
%% for a single channel, concats multiple files into one

% batchf='batch_test'; (name wil be sued to save)
% channel_board=14; % amplifier neural channel (actual baord)

% output: in dir called (CHANNELS_concat), will save metadata, data, and
% copy batch file

%% put names of all files in batch into a cell array

filenames=lt_batchsong_NamesToCell(batchf);


%%  concatenate channels as vectors - save metadata information along the way

fs_all=[];
ampDat_all=[];
metaDat=struct;

% ----- COLLECT DATA
for i=1:length(filenames)
    
    fname=filenames{i};
    
    [amplifier_data,~ ,frequency_parameters, board_adc_data, ...
        board_adc_channels, amplifier_channels, ~] = pj_readIntanNoGui(fname);
    
    % --- collect all amp dat
    ind=[amplifier_channels.chip_channel]==channel_board;
    dattmp=amplifier_data(ind, :); % dat for this chan
    
    ampDat_all=[ampDat_all dattmp]; % concat
    
    % --- collect all fs
    fs_all=[fs_all frequency_parameters.amplifier_sample_rate];
    
    % --- collect metadata
    metaDat(i).filename=fname; % filename
    metaDat(i).numSamps=length(dattmp); % length of file (samps)
    metaDat(i).fs=frequency_parameters.amplifier_sample_rate; % fs
    metaDat(i).songDat=board_adc_data(1,:);
end


% --- confirm that all have same sample rate
if length(fs_all)>1
    assert(all(diff(fs_all)==0), 'problem, FS not same for all file')
end


%% === saving

% --- save a data file with name corresponding
dirsavename=['Chan' num2str(channel_board) 'amp-' batchf];

currdir=pwd;

% --- make a new dir for this data
if ~exist(dirsavename, 'dir')
    % make dir
    mkdir(dirsavename)
    cd(dirsavename)
    
    % --- save
    data=ampDat_all; % wave_clus needs this name change
    save(['data.mat'], 'data');  % save vector
    save(['MetaDat.mat'], 'metaDat'); % save metadata
    eval(['!cp ../' batchf ' .']) % copy batchfile over
else
    if strcmp(input('data already exists - overwrite? (y or n) ' , 's'), 'y')
        eval(['!rm -r ' dirsavename ]);
        mkdir(dirsavename)
        cd(dirsavename)
        
        disp(['... deleted folder ' dirsavename]);
        
        % --- save
        data=ampDat_all; % wave_clus needs this name change
        save(['data.mat'], 'data');  % save vector
        save(['MetaDat.mat'], 'metaDat'); % save metadata
        eval(['!cp ../' batchf ' .']) % copy batchfile over
    else
        % then don't save data
    end
end

cd(currdir)


