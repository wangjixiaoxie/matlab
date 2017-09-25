function nsyls = lv_countsegments(batchname)

fid = fopen(batchname,'r');
[~,numlines] = fscanf(fid,'%s');
frewind(fid);

nsyls = 0;

for i=1:numlines
    fn=fgetl(fid);
    load([fn '.not.mat']);
    nsyls = nsyls+length(labels);
end