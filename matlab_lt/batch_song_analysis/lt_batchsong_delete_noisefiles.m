%% LT 12/4/15 - run in day folder. uses cleandirAuto and deletes all nonsongs
function lt_batchsong_delete_noisefiles(quick_mode, rerun_cleandir)

% quick_mode=1; % if 1, does not prompt user, just goes ahead and deletes
% rerun_cleandir=1; % reruns cleanDirAuto, using more conservative params


%%  check if have run code in this folder before - if has, then promts user to not continue
tmp=dir('DONE_lt_batchsong_delete_noisefiles');
if ~isempty(tmp);
    disp('-- NOTE: this code has been run in this folder before. Might not want to refilter songs');
    DoContinue=input('Continue? (1 or 0) ');
    
assert(DoContinue==1, 'STOPPED');

end
    

%% 

% === find all nonsongs
if rerun_cleandir==0;
    tmp=dir('batch.dcrd');
    
    if isempty(tmp)
        % then have not yet run cleandir
        lt_make_batch(1);
    end
else
    lt_make_batch(3);
    cleandirAuto('batch',1000, 4,4);
    
end

% === delete nonsongs
if quick_mode==1;
    do_delete=1;
else
    % -- first show random subset of files that will be deleted
    lt_disp_random_songs_syls('batch.dcrd',0,1,48)
    edit batch.dcrd
    
    % ask user
    do_delete=input('delete all files from batch.dcrd? (1 or 0) (2 to show more dcrd songs) ');
    
end

if do_delete==2;
    while do_delete==2;
    close all;
     lt_disp_random_songs_syls('batch.dcrd',0,1,48);
   
    % ask user
    do_delete=input('delete all files from batch.dcrd? (1 or 0) (2 to show more dcrd songs) ');
    end
end
   

if do_delete==1;
    
    disp(' --- DELETING FROM batch.dcrd: ');
    
    
       % === write a file indicating that this code was run in this folder
   fid_record=fopen('DONE_lt_batchsong_delete_noisefiles' ,'w');
        fprintf(fid_record,'%s\n','FILES DELETED:');

    % --- open batch and copy all contents to a record file
    fid=fopen('batch.dcrd');    
    fline=fgetl(fid);
    
    while ischar(fline)
       
        % ==== make a record of files deleted
        fprintf(fid_record,'%s\n',fline);

%         % delete this file
%         eval(['!rm ' fline]);
%         eval(['!rm ' fline(1:end-5) '.rec'])
        
        disp(fline)
        
        fline=fgetl(fid);
    end
    
    % ==== actually delete those files
    eval('!rm $(less batch.dcrd)'); % remove all files in that batch file
    
    % -- DONE
   disp('--- DONE');
   tstamp=lt_get_timestamp(0);
        fprintf(fid_record,'%s\n',['PERFORMED ON: ' tstamp]);
        
   fclose(fid_record);
   
else
    disp('DID NOT DELETE ANYTHING!!');
end
