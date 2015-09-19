function [fvals,tvals,pitchcurves]=summary_stats2010(batchnote,filetype,FFwindow,note)
% e.g. summary_stats2010('batchnotes',1,[2000 3000],'a'); 
%       filetype=1 denotes *.cbin
%       fileype=0 denotes *.wav
if filetype==0
    type='wav';
    fs=44100;
else
    type='obs0';
    fs=32000;
end
fvals=findwnoteJC(batchnote,note,'','',0,[2000 2700],8500,1,type,1);
tvals=timing3(fvals);
%shifted=jc_AlignCT(fvals,toff,ind);
for i=1:length(fvals)
    shifted(i,:)=fvals(i).datt;
end
pitchcurves=jc_pitchmat1024(shifted,1024,1020,1,FFwindow(1),FFwindow(2),[1],type,1);
