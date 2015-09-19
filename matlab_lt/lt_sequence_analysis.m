%% 09/19/13
% for pu13bk43, did 7/16, 7/17, label all songs, very conservatively
% and look at labels.  do below, open labels... file, and see that much more
% of the canonical motif is in the WN days, (i.e. replace motif with a single LETTER).


clear all; close all
% batch=input('batch name? ','s');
batch='batch.catch.keep_AllSylLabeled_091913';
fid=fopen(batch);

clock=1;
fid_labels=fopen('labels_compiled','w');
while 1
    fn=fgetl(fid);
    if ~ischar(fn),
    break; 
    end
    temp=open([fn '.not.mat']);
    labels{clock}=temp.labels;
    fprintf(fid_labels,'%s\n',temp.labels);
    clock=clock+1;
%     suffix_start_index=findstr(fn,'.cbin');
%     rec_name_from_cbin=[fn(1:suffix_start_index-1) '.rec'];
%     rec_files_list(clock)=dir(rec_name_from_cbin);
%     clock=clock+1;
end
fclose(fid_labels);


