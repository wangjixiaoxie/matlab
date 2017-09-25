function lv_resegment(batchname,filetype)

fid = fopen(batchname,'r');
[~,numlines] = fscanf(fid,'%s');
frewind(fid);

for i=1:numlines
    fn{i}=fgetl(fid);
end

    
% SegmentFiles(fn,25,400,23,150,filetype); %Minimum + Maximum syl lengths; gap lengtha; gap wqar 3, syl 150;
SegmentFiles(fn,25,150,3,150,filetype); %Minimum + Maximum syl lengths; gap lengtha; gap wqar 3, syl 150;
