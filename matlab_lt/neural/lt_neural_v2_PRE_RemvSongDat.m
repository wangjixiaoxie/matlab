function lt_neural_v2_PRE_RemvSongDat(SummaryStruct)
%% remove song dat from metadat


% THIS IS UP TO DATE. WILL SAVE SONG IN UPPER DIR AND REMOVE SONG DAT FROM
% META DAT. RUNNING WILL NOT CORRUPT ANYTHING. 
% RUN THIS FOR NEW FINALIZED NEURONS.

numbirds = length(SummaryStruct.birds);
for i=1:numbirds
    
    numneurons = length(SummaryStruct.birds(i).neurons);
    
    for ii=1:numneurons
        disp(['========= bird ' num2str(i) ' neuron ' num2str(ii) ', ' ...
            SummaryStruct.birds(i).neurons(ii).batchfilename ' (if blank did nothing)']);
        
        Datstruct = SummaryStruct.birds(i).neurons(ii);
        cd(Datstruct.dirname)
        
        % 1) === IS THIS BATCH'S SONG ALREADY SAVED UP ONE LEVEL?
        cd ..
        songfname = ['SongDat_' Datstruct.batchfilename '.mat'];
        songAlreadyExists = 0;
        if exist(songfname, 'file')==2
            % then already exists, skip
            songAlreadyExists = 1;
        end
        cd(Datstruct.dirname)
        
        % 2) ====  Extract and save song if necessary.
        if songAlreadyExists == 0
            disp('EXTRACTING SONG ...');
            % -- if not saved, then save
            metadat = load('MetaDat');
            
            % save song data in cell array
            numsongs = length(metadat.metaDat);
            SongCellArray = {};
            for j=1:numsongs
                SongCellArray = [SongCellArray single(metadat.metaDat(j).songDat)];
            end
            
            % save
            cd ..
            save(songfname, 'SongCellArray', '-v7.3');
            disp(['-- extracted and saved ! (' songfname ')']);
            cd(Datstruct.dirname)            
        end
        
        % 3) ===== if Song not removed from MetaDat, then do that.
        if exist('DONE_RemovedSongDatFromMetaDat', 'file') ==2
            % then already removed, DO NOTHING.
            
        else
            disp('REMOVING SONG FROM METADAT');
            metadat = load('MetaDat');
            % ------ 2) don't yet remove from metaDat (for backwards compatibility testing)
            % INSTEAD CHANGE THE NAME. IF WORKS FINE, THEN REMOVE.
            metaDat = rmfield(metadat.metaDat, 'songDat'); % also name change
            
            % move old metadat to new file name
            % DELETE OLD METADAT
            eval('!rm MetaDat.mat');
            disp('DELETED OLD METADAT!');
            
            % save new
            save('MetaDat.mat', 'metaDat');
            
            fid = fopen('DONE_RemovedSongDatFromMetaDat', 'w');
            fclose(fid);
        end
    end
end
disp('---------- DONE!');

