function lt_batch_combine(NewBatchSuffix,varargin);

%% 4/19/15 LT - takes multiple song batches and combines and sorts by date and returns one new batch
% e.g. lt_batch_combine('CombinedBatch','batch1','batch2');


%%

for i=1:length(varargin);
    InputBatches{i}=varargin{i};
end

% start new combined batch file

fidnew=fopen(['batch.combined.' NewBatchSuffix],'w');



%%

for i=1:length(InputBatches);
    
    fid=fopen(InputBatches{i},'r');
    fn=fgetl(fid);
    while ischar(fn);
        
        fprintf(fidnew,'%s\n',fn);
        fn=fgetl(fid);
        
    end
    
    fclose(fid);
end

fclose(fidnew);

% Sort by date
lt_sort_batch(['batch.combined.' NewBatchSuffix],'date');






