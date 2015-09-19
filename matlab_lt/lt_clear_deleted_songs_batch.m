function lt_clear_deleted_songs_batch(batch)
%% LT 7/25/15 - removes song filenemas from batch file if the .cbin file is missing from the folder
% INPUTS:
% 1) leave empty to do this for all batches in folder: [lt_clear_deleted_songs_batch('')
% 2) enter batch name to do for only one batch

if isempty(batch)
    do_all_batch=1;
else
    do_all_batch=0;
end


%% ===================
if do_all_batch==0;
    % ========== RUN ON SINGLE DESIGNATED BATCH
    tstamp=lt_get_timestamp(0);
    
    fid=fopen(batch); % old batch
    new_batch=[batch '.' tstamp];
    fidnew=fopen(new_batch, 'w'); % new batch
    
    line=fgetl(fid);
    
    while ischar(line);
        
        % check if line (song file) exists in this folder
        test=dir([line]);
        
        if ~isempty(test);
            % then write this song to new batch file
            fprintf(fidnew, '%s\n', line);
        end
        
        
        line=fgetl(fid);
    end
    
    % close both files and overwrite old batch with new batch(i.e. with batch
    % with deleted songs removed)
    fclose(fid);
    fclose(fidnew);
    
    eval(['!mv '  new_batch ' ' batch]);
    
elseif do_all_batch==1;
     % =================================== DO ALL BATCHES
   
    % ---- COLLECT ALL BATCHES
    batches_all=dir('batch*');
    for k=1:length(batches_all);
        batch=batches_all(k).name;

        % ==== RUN AS ABOVE
        tstamp=lt_get_timestamp(0);
        
        fid=fopen(batch); % old batch
        new_batch=[batch '.' tstamp];
        fidnew=fopen(new_batch, 'w'); % new batch
        
        line=fgetl(fid);
        
        while ischar(line);
            
            % check if line (song file) exists in this folder
            test=dir([line]);
            
            if ~isempty(test);
                % then write this song to new batch file
                fprintf(fidnew, '%s\n', line);
            end
            
            
            line=fgetl(fid);
        end
        
        % close both files and overwrite old batch with new batch(i.e. with batch
        % with deleted songs removed)
        fclose(fid);
        fclose(fidnew);
        
        eval(['!mv '  new_batch ' ' batch]);
    end
end
