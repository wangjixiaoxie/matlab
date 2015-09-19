%% LT 2/13/15 - takes a batch as input, outputs batch.LABELED and batch.UNLABELED
function lt_sort_batch_by_labeled(batch);
%
% e.g.
% batch = 'batch.keep';

fid=fopen(batch,'r'); % open

fn=fgetl(fid); % song name, line by line

fid_l=fopen([batch '.LABELED'],'w');
fid_u=fopen([batch '.UNLABELED'],'w');

while fn~=-1; % -1 means no more lines in file
    
    fn_notmat=[fn '.not.mat']; % go from .cbin to .not.mat file name
    
    if exist(fn_notmat,'file') % does the fn_notmat file exist?
        
        NM=load(fn_notmat);  % load .not.mat file
        
        notes=regexp(NM.labels, '\w', 'once'); % look for first thing that is letter (low or up) or num
        if ~isempty(notes);  % if was labeled
            fprintf(fid_l,'%s\n',fn);
        else % if not labeled
            fprintf(fid_u,'%s\n',fn);
        end
        
    else % if not.mat file does not exist, then obviously cannot be labeled
        fprintf(fid_u,'%s\n',fn);
    end
    
    fn=fgetl(fid);
end


fclose(fid);
fclose(fid_l);
fclose(fid_u);
end
