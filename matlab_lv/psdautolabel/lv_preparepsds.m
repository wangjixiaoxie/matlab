function lv_preparepsds(batchname,filetype,npsds)
%calculates psds for batchfile, should already be segmented/notmat file
%use lv_resegment before

fid = fopen(batchname,'r');
[~,numlines] = fscanf(fid,'%s');
frewind(fid);

for i=1:numlines
    fn{i}=fgetl(fid);
end


[PSDSylOut, SylOut, alllabels, Fs, gapdur]=SegmentSylNormPSD(fn,filetype,npsds);

save(['psddata_' num2str(npsds) 'parts_' batchname '.mat'],'PSDSylOut','SylOut','alllabels','Fs','gapdur')
