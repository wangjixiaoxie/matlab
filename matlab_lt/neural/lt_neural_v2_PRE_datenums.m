function lt_neural_v2_PRE_datenums(SummaryStruct)

%% lt 8/10/17 - put datenums into metadat (since takes a long time) permanent


numbirds = length(SummaryStruct.birds);
for i=1:numbirds
    
    numneurons = length(SummaryStruct.birds(i).neurons);
    
    for ii=1:numneurons
%         disp(['========= bird ' num2str(i) ' neuron ' num2str(ii) ', ' ...
%             SummaryStruct.birds(i).neurons(ii).batchfilename]);
        
        Datstruct = SummaryStruct.birds(i).neurons(ii);
        cd(Datstruct.dirname);
        
        tmp = load('MetaDat.mat');
        if isfield(tmp.metaDat, 'song_datenum')
           continue
        end
        
        % ---- get song dates and put into metadat
%         disp('extracting datenum to metadat');
        numsongs = length(tmp.metaDat);
        for j=1:numsongs
           tmp.metaDat(j).song_datenum = lt_neural_fn2datenum(tmp.metaDat(j).filename);
        end
        
        % -- save
        metaDat = tmp.metaDat;
        save('MetaDat.mat', 'metaDat')
        
    end
end

