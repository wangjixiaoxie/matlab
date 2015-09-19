function [] = db_change_syllable_in_batchfile( batchfile, original_syllable, new_syllable )
%db_change_syllable_in_batchfile Input a batchfile, original_syllable, and
%new_syllable, and it will change all the 'labels' variables in the
%.not.mat files to replace the old syllables with the new.

%opens the batchfile and gets the first line
fid = fopen(batchfile,'r');
current_line = fgetl(fid);

if length(original_syllable) == length(new_syllable)
    motif_length = length(original_syllable)-1;
else
    display('Original syllable and new syllable are different lengths')
    return
end

while ischar(current_line)
    %see if batchfile is a .not.mat file or a .cbin/.wav file, then loads
    %the first line .not.mat file
    if strfind(current_line,'.not.mat')
        load(current_line)
    else
        try
        load([current_line '.not.mat'])
        catch
            current_line = fgetl(fid);
            continue
        end
    end
    
    %searchs in the labels variable for the original syllable and replaces
    %it with the new one (for motif length of 0, it just indexes (much
    %faster))
    if motif_length == 0
        labels(strfind(labels,original_syllable)) = new_syllable;
    else
        matches = strfind(labels,original_syllable);
        for i = 1:length(matches)
            labels(matches(i):matches(i)+motif_length) = new_syllable;
        end
    end
    
    
    %saves the .not.mat file with the modified labels variable
    if strfind(current_line,'.not.mat')
        save(current_line,'Fs','fname','labels','min_dur','min_int','offsets','onsets','sm_win','threshold');
    else
        save([current_line '.not.mat'],'Fs','fname','labels','min_dur','min_int','offsets','onsets','sm_win','threshold');
    end
    
    %goes to the next line in the batchfile
    current_line = fgetl(fid);
end

%close the batch file
fclose(fid);
    


end

