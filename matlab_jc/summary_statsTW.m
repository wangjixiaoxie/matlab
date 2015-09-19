function [shifted,FFstats,avg,toffset,pitch]=summary_statsTW(batch,batchnote)
fs=32000;
fvals=findwnoteJC(batchnote,'b','','',0,[2000 3000],8500,1,'obs0',0);

[vals,trigs]=triglabel(batch,'b',1,1,0,1);
toff=[];
for ii=1:length(trigs)
toff=[toff;trigs(ii).toffset];
end
stdtoff=std(toff);

vvals=getvals(fvals,1,'TRIG');
FFstats=vvals(:,2);


%shifted=jc_AlignCT(fvals,toff,ind);
for i=1:length(fvals)
    shifted(i,:)=fvals(i).datt;
end
[pitch,avg]=jc_pitchmat1024(shifted,1024,1020,2,2000,3000,[1 2],'obs0',1);
toffset=((toff/1000)*(fs)-512)/4+240; %add 160 b/c use findwnoteJC
figure;hold on;
plot(avg);
