dirf('*.cbin','batch')
findcatch('batch')
evsonganaly
edit batch.catch

fvals614Amiss=findwnoteJC('batch614Anotes','b','','',0,[2000 2700],8500,1,'obs0',1);
for i=1:length(fvals614Amiss)
    shifted614Amiss(i,:)=fvals614Amiss(i).datt;
end
pitch614Amiss=jc_pitchmat1024(shifted614Amiss,1024,1020,2,2000,2700,[1],'obs0',1);
fvals614Ahit=findwnoteJC('batch614Anotes','a','','',0,[2000 2700],8500,1,'obs0',1);
for i=1:length(fvals614Ahit)
    shifted614Ahit(i,:)=fvals614Ahit(i).datt;
end
pitch614Ahit=jc_pitchmat1024(shifted614Ahit,1024,1020,2,2000,2700,[1],'obs0',1);
pitch614A=[pitch614Amiss,pitch614Ahit];
