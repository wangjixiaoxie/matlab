function [  ] = db_change_miss_label(batchfile, motif, new_syllable, syl_position )
%db_change_miss_label Relabels syllable to new_syllable in a batchfile if
%the syllable was a miss (so you can make a new template to misses)
%   db_change_miss_label(batchfile, syllable, new_syllable)
%   Loads the .rec and .not.mat files, looks for the onset and offset times
%   of labeled syllables. Then sees if the hit times in the .rec file fall
%   between those times. If so, syllable stays the same, else it changes it
%   to new_syllable.

%% Prepares motif/syllable search pattern
if nargin < 4
    syl_position = 0;
end


%% Changes miss labels
%opens the batchfile
fid_batch = fopen(batchfile,'r');
tline = fgetl(fid_batch);

while ischar(tline)
    try
        %loads the .not.mat file
        if ~isempty(strfind('.not.mat',tline))
            load tline
            tline = tline(1:end-8);
        else
            load([tline '.not.mat'])
        end
        
        %reads the .rec file
        fid_rec = fopen([tline(1:end-4) 'rec'],'r');
        j = 1;
        hit_time = {};
        nextline = fgetl(fid_rec);
        %goes through .rec file and records hit times
        while ischar(nextline)
            if ~isempty(strfind(nextline, 'msec:'))
                hit_time{j} = str2double(nextline(1:6))*1000;
                j = j+1;
            end
            nextline = fgetl(fid_rec);
        end
        fclose(fid_rec);
    
        %finds syllable positions of interest
        positions = regexp(labels, motif) + syl_position;
        %finds time when syllables occured
        time_of_syl = [onsets(positions) offsets(positions)];
        
        % if motif is not found, then song is skipped
        if isempty(positions)
            display([tline ' did not have syllable of interest'])
        % if motif is found but no hits, then all labels are changed
        elseif ~isempty(positions) && isempty(hit_time)
            labels(positions) = new_syllable;
        else
            hit_time = cell2mat(hit_time)';
            for i = 1:size(time_of_syl,1)
                if sum(time_of_syl(i,1) < hit_time & time_of_syl(i,2) > hit_time) == 0
                    labels(positions(i)) = new_syllable;
                end
            end
        end
                
        save([tline '.not.mat'], 'labels', '-append')
        display(tline)
    catch err
        display([tline ' has not .not.mat file'])
    end
    tline = fgetl(fid_batch);
end
fclose(fid_batch);


end

