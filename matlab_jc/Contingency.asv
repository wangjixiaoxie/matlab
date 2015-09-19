function Contingency(toff,batchNOTE)
fvals=findwnote4(batchNOTE,'a','','',0,[2000 2700],8500,1,'obs0',1);
align_shift=jc_align_shift(fvals);
alignMS=align_shift/32;
figure;hist(toff'-alignMS);
std(toff'-alignMS)