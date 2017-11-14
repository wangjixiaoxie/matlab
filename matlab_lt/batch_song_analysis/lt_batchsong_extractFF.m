function DATSTRUCT = lt_batchsong_extractFF(ListOfDirs_UNDIR, ListOfDirs_DIR, ListOfBatch, MotifsToExtract)


%% lt 11/14/17 - extracts results from lt_batchsong_calcFF

%%
ListOfDirs_ALL = [ListOfDirs_UNDIR ListOfDirs_DIR];

%%
DATSTRUCT = struct;
for mm = 1:length(MotifsToExtract)
    
    regexpr_str = MotifsToExtract{mm};
    DATSTRUCT.motif(mm).motif = regexpr_str;
    
    count=1;
    for i=1:length(ListOfBatch)
        
        dirname = ListOfDirs_ALL{i};
        batchf = ListOfBatch{i};
        
        % ==== is this DIR or UNDIR song?
        if any(ismember(ListOfDirs_UNDIR, dirname))
            isDIR = 0;
        elseif any(ismember(ListOfDirs_DIR, dirname))
            isDIR = 1;
        else
            disp('PROBLEM!!! - dir or undir?');
        end
        
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
            
            % ---------------- extract time of song
            if any(strfind(fname, '.rhd'));
                [dtnum datestring]=lt_neural_fn2datenum(fname);
            else
                disp('PROBLEM!!! - name extraction, do for evtaf');
                asdfsdf;
            end
            
            % =============== find all instances of this motif in this
            % songfile
            [tokenExtents, startinds, endinds, matchlabs] = lt_batchsong_regexp(notmat.labels, regexpr_str);           %%
            if isempty(tokenExtents)
                fname = fgetl(fid);
                continue
            end
            assert(size(tokenExtents,1)==1, 'needs to be horizontal');
            
            for j=tokenExtents
               
                % ============= extract data for this rend
                    DATSTRUCT.motif(mm).rendnum(count).ff = calcff.FFstruct.FFall(j);
                    DATSTRUCT.motif(mm).rendnum(count).syl = notmat.labels(j);
                    DATSTRUCT.motif(mm).rendnum(count).fname = fname;
                    DATSTRUCT.motif(mm).rendnum(count).dirname = dirname;
                    DATSTRUCT.motif(mm).rendnum(count).datenum_song_SecRes = dtnum;
                    DATSTRUCT.motif(mm).rendnum(count).datestr = datestring;
                    DATSTRUCT.motif(mm).rendnum(count).time_withinsong = notmat.onsets(j)/1000;
                    DATSTRUCT.motif(mm).rendnum(count).isDIR = isDIR;
                    count = count+1;
            end
            
            
            fname = fgetl(fid);
        end
    end
    
end
