function [fvals,FFstats,toffset]=summary_statsTW2(batch,batchnote)
fs=32000;
fvals=findwnoteJC(batchnote,'l','','',0,[3200 3800],8500,1,'obs0',0);

[vals,trigs]=triglabel(batch,'l',1,1,0,0);
toff=[];
for ii=1:length(trigs)
toff=[toff;trigs(ii).toffset];
end
stdtoff=std(toff);

vvals=getvals(fvals,1,'TRIG');
FFstats=vvals(:,2);

%toffset=((toff/1000)*(fs)-512)/4; %add 160 b/c use findwnoteJC
%shifted=jc_AlignCT(fvals,toff,ind);

toffset=((toff/1000)*(fs)-512)/4+240;


