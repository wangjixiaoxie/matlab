function [outdir] = lt_neural_SnapshotCurrDat(SummaryStruct)
%% lt 9/12/17 - copies .not.mat to preserve current dataset and labeling.

cd('/bluejay5/lucas/analyses/neural/Snapshot_Dat_Labels');

tstamp = lt_get_timestamp(0);
savedir = tstamp;

% 1) save summarystruct
mkdir(savedir);
cd(savedir);
save('SummaryStruct', 'SummaryStruct');

numbirds = length(SummaryStruct.birds);

for i=1:numbirds
    
    mkdir(['bird' num2str(i)]);
    cd(['bird' num2str(i)]);
    numneurons = length(SummaryStruct.birds(i).neurons);
    
    for ii=1:numneurons
        
        mkdir(['neuron' num2str(ii)]);
        cd(['neuron' num2str(ii)]);
        
        % ---- will move stuff to this dir
        currdir = pwd;
        
        % ---- copy data over
       batchf = SummaryStruct.birds(i).neurons(ii).batchfilename;
       dirname = SummaryStruct.birds(i).neurons(ii).dirname;
       
       cd(dirname);
       
       % ==== save FF extracted.
       if exist('extractFF_params.mat')==2
       eval(['!cp extractFF_params.mat ' currdir]);
       eval(['!cp extractFF.mat ' currdir]);
       end
       
       
       % ===== for each file in batch file copy that batch file to a save
       % area
       cd ..
       fid = fopen(batchf);
       fline = fgetl(fid);
       while ischar(fline)
           
           notmatfname = [fline '.not.mat'];
           
           if exist(notmatfname)~=2
               disp(['skipping ' notmatfname ' (not exist)']);
           else
           eval(['!cp ' notmatfname ' ' currdir]);
           disp(['copying ' notmatfname]);
           end
           fline = fgetl(fid);
       end
       
       cd(currdir);
       cd ..
        
    end
    
    cd ..
end

outdir = ['/bluejay5/lucas/analyses/neural/Snapshot_Dat_Labels/' tstamp]
