[vals,trigs]=triglabel('batchJCfiles','a',1,1,0,0);
toff=[];
for ii=1:length(trigs)
toff=[toff;trigs(ii).toffset];
end
stdtoff=std(toff);
fs=32000;
toffset=240+((toff/1000)*(fs)-512)/4;
fvalsampon=findwnoteJC('batchJCnotes','c','','',0,[2000 2700],8500,1,'obs0',1);
for i=1:length(fvalsampon)
    shiftedampon(i,:)=fvalsampon(i).datt;
end
pitchampon=jc_pitchmat1024(shiftedampon,1024,1020,2,2000,2700,[1],'obs0',1);