figure
ps.ax=subplot(211);
ps.minx=5
ps.maxx=5
ps.ac_col='k'
ps.mu_col=[0.4 0.4 1]
ps.plotarrow=0
ps.flip=1
ps.plotraw=1
ps.axraw=subplot(212);
ps.addzero=0
ps.plotsep=1;
ps.type='nor';
ps.insert=1
ps.aligntimes=1
ps.subdayinterp=1
ps.usepct=1
figure

[outvlaczcomb,outvlmuzcomb,sdind]=plotasympdynamics(sumdynasymp,ps)
