%% lt 7/7/17 - prints labels for all songs in a batch
function lt_batch_disp_all_labels(batchname)

fid = fopen(batchname);

fline = fgetl(fid);

AllSyls = [];
while ischar(fline)
    
%     dots = strfind(fline, '.');
    
    fname_notmat = [fline '.not.mat'];
    
    if exist(fname_notmat, 'file') ==2;
        tmp = load(fname_notmat);
        disp([fname_notmat ': ' tmp.labels]);
        AllSyls = [AllSyls tmp.labels];
    else
        disp([fname_notmat ': not labeled']);
    end
    
   
    fline = fgetl(fid);
end

sylsused = unique(AllSyls);
disp(' ');
disp(['SYLS USED: ' sylsused]);