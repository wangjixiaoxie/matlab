%short pulse data list.m
clear sps
%this is for bk48w74
sls(1).inds=[3]

sls(1).stlen=[60 ] 
sls(1).datbnds{1}=[0 0.09]
sls(1).fb{1}=1
sls(1).catch{1}=1;
sls(1).OFF_TMS{1}=[.07 .078]
%  sps(1).OFF_TMS{2}=[.075 .083]
% sps(1).OFF_TMS{3}=[.081 .089]
% sps(1).OFF_TMS{3}=[.081 .09]
sls(1).RESIDTMS=[.07 .1];
sls(1).stind=1
sls(1).drxn='do'

