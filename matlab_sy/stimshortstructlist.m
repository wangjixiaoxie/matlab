%short pulse data list.m
clear sps
%this is for bk48w74
sps(1).inds=[4 8 9]
sps(1).datbnds{1}=[0 0.016]
sps(1).datbnds{2}=[0 .016]
sps(1).datbnds{3}=[.03 .046]

sps(1).stlen=[60 5 5] 
sps(1).catch{1}=[1 2 3]
sps(1).catch{2}=[1 2 3]
sps(1).catch{3}=[1 2 3]
sps(1).fb{1}=1
sps(1).fb{2}=2
sps(1).fb{3}=3
sps(1).OFF_TMS{1}=[.074 .082]
 sps(1).OFF_TMS{2}=[.076 .084]
sps(1).OFF_TMS{3}=[.079 .087]
% sps(1).OFF_TMS{3}=[.081 .09]
sps(1).RESIDTMS=[.07 .1];
sps(1).stind=1
sps(1).drxn='do'

sps(2).inds=[30] 
% sps(2).datbnds{1}=[.003 .019]
sps(2).datbnds{1}=[.05 .066]
sps(2).stlen=[ 10]
sps(2).catch{1}=[1 ]
sps(2).catch{2}=[1 ]
sps(2).fb{1}=1
% sps(2).fb{2}=2
sps(2).stind=1
sps(2).OFF_TMS=sps(1).OFF_TMS
sps(2).RESIDTMS=[.07 .105]
sps(2).drxn='up'
%bk59w37

% sps(3).inds=[29 30 10 11]
sps(3).inds=[29 30]
sps(3).datbnds{1}=[.028 .041]
sps(3).datbnds{2}=[.041 .054]
sps(3).datbnds{3}=[.028 .054]
% sps(3).datbnds{3}=[.012 .028]
% sps(3).datbnds{4}=[.012 .028]
sps(3).stlen=[20 20]
sps(3).catch{1}=[1 2 ]
sps(3).catch{2}=[1 2 ]
sps(3).catch{3}=[1 2 ]
% sps(3).catch{3}=[1 2 ]
% sps(3).catch{4}=4
sps(3).fb{1}=[1 2]
sps(3).fb{2}=[1 2]
sps(3).fb{3}=[1 2]
sps(3).drxn='up'

% sps(3).fb{3}=3
% sps(3).fb{4}=4              
sps(3).stind=2
sps(3).OFF_TMS{1}=[.073 .081]
% sps(3).OFF_TMS{2}=[.078 .088]
% sps(3).OFF_TMS{3}=[.083 .093]
sps(3).RESIDTMS=[.075 .105];



sps(4).inds=[31 21]
sps(4).datbnds{1}=[.023 .036]
sps(4).datbnds{2}=[.014 .030]
sps(4).stlen=[5 40] 
sps(4).catch{1}=[1 2]
sps(4).catch{2}=[1 2]
sps(4).catch{3}=[1 2]
sps(4).fb{1}=1
sps(4).fb{2}=2
sps(4).drxn='do'
sps(4).stind=2
sps(4).OFF_TMS=sps(3).OFF_TMS;
sps(4).RESIDTMS=sps(3).RESIDTMS;



% %bk57w35orig.              
% sps(5).inds=[45 35 36 [46:51]  41]
% sps(5).stlen=[15 15 40 15 15 60]
% sps(5).datbnds{1}=[.054 .070]
% sps(5).datbnds{2}=[.044 .060]
% sps(5).datbnds{3}=[.027 .043]
% sps(5).datbnds{4}=[.048 .063]
% sps(5).datbnds{5}=[.063 .079]
% sps(5).datbnds{6}=[.03 .046]
% sps(5).catch{1}=[1 2 3]
% sps(5).catch{2}=[1 2 3]
% sps(5).catch{3}=[1 2 3]
% sps(5).catch{4}=[4:10]
% sps(5).catch{5}=[4:10]
% sps(5).catch{6}=[4:10]
% sps(5).fb{1}=1
% sps(5).fb{2}=2
% sps(5).fb{3}=3
% sps(5).fb{4}=4:9
% sps(5).fb{5}=4:9
% sps(5).fb{6}=10
% sps(5).stind=3
% sps(5).OFF_TMS{1}=[.075 .085]
% sps(5).OFF_TMS{2}=[.08 .09]
% sps(5).OFF_TMS{3}=[.085 .095]
% sps(5).RESIDTMS=[.075 .12]

%bk57w35orig.              
sps(5).inds=[45 35  46:51 ]
sps(5).stlen=[15 15 15 15 15 15 15 15 ]
sps(5).datbnds{1}=[.054 .070]
sps(5).datbnds{2}=[.044 .060]
sps(5).datbnds{3}=[.048 .063]
sps(5).datbnds{4}=[.063 .079]
sps(5).datbnds{5}=[.044 .079]
% sps(5).datbnds{5}=[.03 .046]
sps(5).catch{1}=[1 2 ]
sps(5).catch{2}=[1 2 ]
sps(5).catch{3}=[3:8]
sps(5).catch{4}=[3:8]
sps(5).catch{5}=[1:8]
% sps(5).catch{5}=[3:8]

sps(5).fb{1}=1
sps(5).fb{2}=2
sps(5).fb{3}=3:8
sps(5).fb{4}=3:8
sps(5).fb{5}=1:8

sps(5).stind=3
sps(5).OFF_TMS{1}=[.075 .085]
sps(5).OFF_TMS{2}=[.08 .09]
sps(5).OFF_TMS{3}=[.085 .095]
sps(5).RESIDTMS=[.075 .12]
sps(5).drxn='do'





%  %bk57w35evren 
% sps(6).inds=[35 36]
% sps(6).stlen=[10 10 60 60]
% sps(6).datbnds{1}=[.026 .042]
% sps(6).datbnds{2}=[.042 .056]
% sps(6).datbnds{3}=[.007 .023] 
% sps(6).datbnds{4}=[.019 .035]
% sps(6).catch{1}=[1 2]
% sps(6).catch{2}=[1 2]
% sps(6).catch{3}=[1 2]
% sps(6).catch{4}=[1 2]
% sps(6).fb{1}=[1]
% sps(6).fb{2}=[1]
% sps(6).fb{3}=2
% sps(6).fb{4}=2
% sps(6).stind=4
% sps(6).OFF_TMS=sps(5).OFF_TMS;
% sps(6).RESIDTMS=sps(5).RESIDTMS;

 %bk57w35evren 
sps(6).inds=[35 ]
sps(6).stlen=[10  ]
sps(6).datbnds{1}=[.026 .042]
sps(6).datbnds{2}=[.042 .056]
sps(6).datbnds{3}=[.026 .056]
% sps(6).datbnds{3}=[.007 .023] 
% sps(6).datbnds{4}=[.019 .035]
sps(6).catch{1}=[1 ]
sps(6).catch{2}=[1 ]
sps(6).catch{3}=[1]
% sps(6).catch{3}=[1 2]
% sps(6).catch{4}=[1 ]
sps(6).fb{1}=[1]
sps(6).fb{2}=[1]
sps(6).fb{3}=[1]
% sps(6).fb{3}=2
% sps(6).fb{4}=2
sps(6).stind=4
sps(6).OFF_TMS=sps(5).OFF_TMS;
sps(6).RESIDTMS=sps(5).RESIDTMS;
sps(6).drxn='up'

sps(7).inds=[10 11 12 ]
sps(7).stlen=[10 10 10 ]
sps(7).datbnds{1}=[.028 .042]
sps(7).datbnds{2}=[.042 .056]
sps(7).datbnds{3}=[.056 .068] 
sps(7).datbnds{4}=[.028 .068]

sps(7).catch{1}=[1 2 3 ]
sps(7).catch{2}=[1 2 3 ]
sps(7).catch{3}=[1 2 3 ]
sps(7).catch{4}=[1 2 3]

sps(7).fb{1}=[1 2 3 ]
sps(7).fb{2}=[1 2 3 ]
sps(7).fb{3}=[1 2 3 ]
sps(7).fb{4}=[1 2 3 ]
sps(7).drxn='do'
sps(7).stind=1
sps(7).OFF_TMS=sps(1).OFF_TMS;
sps(7).RESIDTMS=sps(1).RESIDTMS;

sps(8).inds=[45 46 47 ]
sps(8).stlen=[60 60 10 ]
sps(8).datbnds{1}=[.028 .044]
sps(8).datbnds{2}=[.028 .044]
sps(8).datbnds{3}=[.048 .064] 

sps(8).catch{1}=[1 2 3 ]
sps(8).catch{2}=[1 2 3 ]
sps(8).catch{3}=[1 2 3 ]

sps(8).fb{1}=[1  ]
sps(8).fb{2}=[ 2 ]
sps(8).fb{3}=[  3 ]
sps(8).stind=1
%%%FIX
sps(8).drxn='do'

sps(8).OFF_TMS=sps(1).OFF_TMS;
sps(8).RESIDTMS=sps(1).RESIDTMS;


sps(9).inds=[48 49 50 51 ]
sps(9).stlen=[10  10 10 10 ]
sps(9).datbnds{1}=[.008 .018]
sps(9).datbnds{2}=[0.018 0.028]
sps(9).datbnds{3}=[.028 .038] 
sps(9).datbnds{4}=[.038 .048] 
sps(9).datbnds{5}=[.048 .058] 
sps(9).datbnds{6}=[.058 .068] 
sps(9).datbnds{7}=[.068 .078] 
sps(9).datbnds{8}=[.078 .088]
sps(9).datbnds{9}=[.008 .088]
% sps(9).datbnds{9}=[.071 .079]
% sps(9).datbnds{10}=[.076 .084]
% sps(9).datbnds{11}=[.074 .081]
% sps(9).datbnds{12}=[.0595 .0645]
% sps(9).datbnds{13}=[.0645 .0695]
% sps(9).datbnds{14}=[.0695 .0745]
% sps(9).datbnds{15}=[.0745 .0795]
% sps(9).datbnds{12}=[.082 .089]
% sps(9).datbnds{13}=[.065 .070]
% sps(9).datbnds{14}=[.070 .075]
% sps(9).datbnds{15}=[.080 .085]
% sps(9).datbnds{16}=[.005 .090]




% sps(9).datbnds{10}=[0 0.094];
for ii=1:length(sps(9).datbnds)
sps(9).catch{ii}=[1:4]
end

for ii=1:length(sps(9).datbnds)
sps(9).fb{ii}=[1:4]
end

sps(9).stind=1
sps(9).OFF_TMS=sps(1).OFF_TMS;
sps(9).RESIDTMS=sps(1).RESIDTMS;
sps(9).drxn='do'


sps(10).inds=[65]
sps(10).stlen=[10]
sps(10).datbnds{1}=[.000 .085]
sps(10).fb{1}=1;
sps(10).stind=1;
sps(10).catch{1}=1;
sps(10).OFF_TMS=sps(1).OFF_TMS;
sps(10).RESIDTMS=sps(1).RESIDTMS;

%%%%%%%%%FIX
sps(10).drxn='do'



sps(11).inds=[66]
sps(11).stlen=[10]
sps(11).datbnds{1}=[.000 .085]
sps(11).fb{1}=1;
sps(11).stind=1;
sps(11).catch{1}=1;
sps(11).OFF_TMS=sps(1).OFF_TMS;
sps(11).RESIDTMS=sps(1).RESIDTMS;
%%%%%%FIXX
sps(11).drxn='do'


sps(12).inds=[71 72 73 74]
sps(12).stlen=[10 10 10 10]
sps(12).datbnds{1}=[.000 .015]
sps(12).datbnds{2}=[.010 .025]
sps(12).datbnds{3}=[.020 .035] 
sps(12).datbnds{4}=[.030 .045] 
sps(12).datbnds{5}=[.040 .055] 
sps(12).datbnds{6}=[.050 .065] 
sps(12).datbnds{7}=[.060 .075] 
sps(12).datbnds{8}=[.070 .085] 
sps(12).datbnds{9}=[0 0.08]
for ii=1:9
    sps(12).fb{ii}=1:4;
    sps(12).stind=1;
    sps(12).catch{ii}=1:4;
end
sps(12).OFF_TMS=sps(1).OFF_TMS;
sps(12).RESIDTMS=sps(1).RESIDTMS;
sps(12).drxn='up'
