function [ spec time freq filename ] = lt_db_spec_syllable(batch_file,syllables, pre_time,post_time, pre_note, post_note)

% LT MODIFIED 12/12 to use lt_db_get_avn instead of db_get_avn (line 27.
% difference is that now syl is defined as a fixed distance before offset
% (as opossed to relative to onset). 

%db_spec_syllable creates a structure with the power spectrums of the
%syllables in your batch file
%   


% uses the function db_get_avn to get the power spectrum of every
% syllable in every file. spec_all{i}{j}, where i is
% the file #, and j is the syllable # within that file.
% time are the times of the power spectrums.
% freq are a list of the frequencies

if exist('pre_note','var') == 0
    pre_note = '';
else
end

if exist('post_note','var') == 0
    post_note = '';
else
end


% [spec_all time freq filename] = db_get_avn(batch_file,syllables,pre_time,post_time,pre_note,post_note,'obs0');
[spec_all time freq filename] = lt_db_get_avn(batch_file,syllables,pre_time,post_time,pre_note,post_note,'obs0');

%gets rid of files with no labeled syllables
spec_all = spec_all(~cellfun('isempty', spec_all));
filename = filename(~cellfun('isempty',filename))';

%concatenating filenames into one giant cell
for i = 1:length(filename)
    filename{i} = filename{i}';
end
filename = cat(1,filename{:});

%reorganizes structure so that syllables are all listed in order (no
%separation by file). Will be a 3-d matrix (freq, time, syllable #)
spec = {};
for j = 1:length(spec_all)
    if j == 1
        spec = spec_all{j};
    else
        spec = cat(3,spec,spec_all{j});
    end
end

clear spec_all;

%gets rid of all frequencies above 10 kHz (they are filtered out)
spec = spec(freq<=10000,:,:);
freq = freq(freq<=10000);

%convert time from sec to msec
time = 1000*time;
        



end

