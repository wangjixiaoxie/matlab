
ps.minx=6
ps.maxx=6
ps.ac_col='k'
ps.mu_col=[0.4 0.4 1]

ps.flip=1
ps.plotraw=1
ps.addzero=0
ps.plotsep=1;
ps.type='nor';
ps.insert=1
ps.aligntimes=1
ps.usepct=1
ps.subdayinterp=0
figure

[outvlaczcomb,outvlmuzcomb,sdind]=plotasympdynamics(sumdynasymp(3),ps)
