
ps.minx=5
ps.maxx=5
ps.ac_col='k'
ps.mu_col=[0.4 0.4 1]

ps.flip=1
ps.plotraw=1
ps.addzero=1
ps.plotsep=0;
ps.type='nor';
ps.insert=1
ps.aligntimes=1
ps.usepct=1
figure

[outvlaczcomb,outvlmuzcomb]=plotinitdynamics(sumdyn(1:16),sumbs(1:9),ps)
