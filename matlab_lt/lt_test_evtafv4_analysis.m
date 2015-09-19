clear all; close all;

filename='/newhome/lucast4/Dropbox/SCIENCE/code/MATLAB_notpath/Randy_evtaf4/randy2notethreshold.evconfig2';

config=ReadEvTAFv4ConfigFile(filename);



%% testing on one day
clear all; close all;

cd /bluejay3/lucas/birds/pu37wh20/112414_SeqDepPitchShift2_durWN_day2
batch = 'batch.labeled.all';
config='config_112414_thr.evconfig2';

EvTAFv4Sim_LT(batch,config,'obs0');



%% other useful programs

[ND,OP]=ReadEvTAFv4ConfigFile('config_112414_thr.evconfig2');

WriteEvTAFv4ConfigFile

UIEvTAFv4Sim('batch.labeled.all',ND)

uievtaf4sim('batch.labeled.all',config)
