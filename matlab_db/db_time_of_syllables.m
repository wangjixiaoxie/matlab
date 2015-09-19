function [ motif_onsets, motif_offsets ] = db_time_of_syllables( filename, motif, fs_or_ms )
%db_time_of_syllables Goes through a cbin file and gives the onset and
%offset times for the motif specified. Asks if you want time in
%milliseconds or 1/(sampling frequency)

if isempty(strfind(filename, '.not.mat')) == 1
    load([filename '.not.mat'])
else
    load(filename)
end

motif_length = length(motif)-1;

if strcmpi(fs_or_ms, 'ms')
    motif_onsets = onsets(strfind(labels,motif));
    motif_offsets = offsets(strfind(labels,motif)+motif_length);
elseif strcmpi(fs_or_ms, 'fs')
    motif_onsets = fix(onsets(strfind(labels,motif))*Fs/1000);
    motif_offsets = fix(offsets(strfind(labels,motif)+motif_length)*Fs/1000);
end

end

