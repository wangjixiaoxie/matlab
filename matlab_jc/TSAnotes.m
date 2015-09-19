%Tuesday 10/14/08
%TSA notcatch analysis - analyzed batch13pm (650 notes)
%TSA1014_3.m (does it all)
%sim.m (determines significance)
%saved as TSA_1014.mat in bk50w18 folder
dirf('*.cbin','batch')
evsonganaly
mk_tempf('batch1800',templa9,2,'obs0');
get_trigt2('batch1800',cntrng9,0.3,128,1,1);
vals1800=jctaf_freq('batch1800',[2200 2750],'a',128,'obs0',1,0);
%Appears that lots of recovery has occurred
%AllTSAdata.mat in bk50w18 folder has vals around 10am,1pm(1300)and 6pm(1800)
