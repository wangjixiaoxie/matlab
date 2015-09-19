function [ syl_amp, absolute_time, ratio_amp ] = db_motif_amplitude(cbin_not_mat_file, motif, base_syl )
%db_motif_amplitude Will find the amplitude for each syllable and the amplitude ratio with a base syllable. Input file,
%motif, and the base_syl (if '', will just give amplitude).
%   Will go through each cbin or cbin.not.mat file find the motif of interest. Will
%   then find the amplitude of each syllable (compared to baseline before the song)
%   and the ratio with that and some base syllable (in the motif).

%% Preparing the variables

%Calculates the number of songs to make a cell for timing
[s,song_number] = system(['wc -l ' cbin_not_mat_file]);
clear s
song_number = str2double(song_number(1:strfind(song_number, cbin_not_mat_file)-2));

%Makes a cell that will have when the motif was sung
absolute_time = cell(song_number,1);
syl_amp = cell(song_number,1);
ratio_amp = cell(song_number,1);

%Calculates length of motif
motif_length = length(motif);

%Finds the position of the base_syl
if strcmpi(base_syl,'')
else
    base_syl_location = strfind(motif,base_syl);
end


%% Calculating the duration of the motif
% Opening the cbin_not_mat file
fid = fopen(cbin_not_mat_file, 'r');

%Gets the first line of the input file
next_line = fgetl(fid);
waveform = [];
i = 1;

while ischar(next_line)
    %Checks to see if file is cbin or cbin.not.mat and loads the
    %cbin.not.mat file
    if isempty(strfind(next_line, '.not.mat'))
        load([next_line '.not.mat']);
        raw_waveform = ReadCbinFile(next_line);
    elseif ~isempty(strfind(next_line, '.not.mat'))
        load(next_line)
        raw_waveform = ReadCbinFile(next_line(1:end-8));
    end
    
    %finds the sampling frequency
    sampling = Fs;
    
    %Bandpass filtering the raw_waveform
    low_filter = 500;
    high_filter = 10000;
    
    waveform = bandpass(raw_waveform,Fs,low_filter,high_filter);
    
    
    %finds the rms of the baseline (lowest rms of length 250ms within the
    %first 2 seconds of the cbin file)
    sliding_baseline_window = linspace(0,sampling*2,101);
    baseline_window_length = sampling*.25;
    base_amp = zeros(length(sliding_baseline_window)-1,1);
    for j = 1:length(sliding_baseline_window)-1
        base_amp(j) = rms(waveform(sliding_baseline_window(j)+1:sliding_baseline_window(j)+baseline_window_length));
    end
    base_amp = min(base_amp);
    
    
    for j = 1:motif_length
        try
            %finds the timing of the motif
            occurence_motif(:,j) = [strfind(labels, motif)+j-1]';
            
            absolute_time{i}(:,j) = db_timing4_file(next_line) + onsets(occurence_motif(:,j)).*(1/1000).*(1/86400);
            
            %finds the onset and offset of each syllable in the motif (in
            %sampling time)
            syl_onsets(:,j) = round(onsets(occurence_motif(:,j))*(Fs/1000));
            syl_offsets(:,j) = round(offsets(occurence_motif(:,j))*(Fs/1000));
            for k = 1:size(syl_onsets,1)
                syl_rms(k,j) = rms(waveform(syl_onsets(k,j):syl_offsets(k,j)));
                syl_amp{i}(k,j) = 20*log10(syl_rms(k,j)./base_amp);
                %if in the future want to look amp using 'universal'
                %background noise
                %                 syl_amp(k,j) = 20*log(syl_rms(k,j)./20);
            end
        catch err
            %because the motif did not occur in this song
            syl_amp{i}(1,j) = NaN;
            absolute_time{i}(1,j) = NaN;
        end
    end
    
    if strcmpi(base_syl,'')
    else
        try
            for j = 1:size(syl_amp{i},2)
                ratio_amp{i}(:,j) = syl_amp{i}(:,j)./syl_amp{i}(:,base_syl_location);
            end
        catch err
            ratio_amp{i}(1,j) = NaN;
        end
    end
        
    
    i = i+1;
    next_line = fgetl(fid);
    clear occurence_motif syl_onsets syl_offsets syl_rms
end

fclose(fid);

%Converts motif_timing to a matrix
syl_amp = cell2mat(syl_amp);
ratio_amp = cell2mat(ratio_amp);
absolute_time = cell2mat(absolute_time);

%Gets rid of NaN
syl_amp = syl_amp(~isnan(syl_amp));
ratio_amp = ratio_amp(~isnan(ratio_amp));
absolute_time = absolute_time(~isnan(absolute_time));

%reshapes matrix after isnan function
syl_amp = reshape(syl_amp,length(syl_amp)./motif_length,motif_length);
ratio_amp = reshape(ratio_amp,length(ratio_amp)./motif_length,motif_length);
absolute_time = reshape(absolute_time,length(absolute_time)./motif_length,motif_length);








end

