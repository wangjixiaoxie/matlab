function labelstruct=loadbatchlabels(batchfile);
%labelstruct=loadbatchlabels(batchfile);
%
%labelstruct.fn - filename
%labelsstruct.labels - labels
%

fid=fopen(batchfile,'r');
if (fid==-1)
    error(['batch file ',batchfile,' does not exist']);
end

cnt=0;
while (1)
    fn=fgetl(fid);
    if (~ischar(fn))
        break;
    end
    
    pp=findstr(fn,'.not.mat');
    if (length(pp)==0)
        fn=[fn,'.not.mat'];
    end
    
    if (~exist(fn,'file'))
        continue;
    end
    cnt=cnt+1;
    load(fn);
    labelstruct(cnt).fn=fn;
    labelstruct(cnt).labels=labels;
end
fclose(fid);
return;