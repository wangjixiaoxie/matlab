%make birdstruct
% 
clear sts
sts(1).path='/oriole/bk48w74/datasum'
sts(1).matfilename='datsum.mat'
sts(1).color=[1 .8 .8]
sts(1).basruns=[1 2 19 23 45 65 66 67 ]
sts(1).ntind=1;
sts(1).contrind=[];
sts(1).sfact=[3000 1000];
sts(1).bname='bk48w74'
sts(1).initshiftind=[4 27]
sts(1).skipasymp=0;
sts(1).revasympvl=[];

% bs(1).exsong='figs/
%very difficult to find a control note in this song.

sts(2).path='/oriole6/bk59w37redo/datasum'
sts(2).matfilename='datsum.mat'
sts(2).color='k'
sts(2).ntind=1;
sts(2).contrind=[] ;
sts(2).basruns=[1 17 18];
% sts(1).revruns{1}=[24 27 30 32]
% sts(1).revruns{2}=[46 52 62 65 70]
sts(2).sfact=[3000 2000];
sts(2).bname='bk59w37'
sts(2).initshiftind=[4 19]
sts(2).skipasymp=0;

sts(3).path='/oriole6/bk57w35/datsum'
sts(3).matfilename='datsum.mat'
sts(3).color='k'
sts(3).ntind=1;
sts(3).contrind=[] ;
sts(3).basruns=[ 13 30];
% sts(1).revruns{1}=[24 27 30 32]
% sts(1).revruns{2}=[46 52 62 65 70]
sts(3).sfact=[3000 2000];
sts(3).bname='bk57w35'
sts(3).initshiftind=[3 22 34]
sts(3).skipasymp=0;

sts(4).path='/oriole2/bk57w35evren/bk57w35/datasum'
sts(4).matfilename='datsum.mat'
sts(4).color='k'
sts(4).ntind=1;
sts(4).contrind=[] ;
sts(4).basruns=[13 14 28 46]; %47];
% sts(1).revruns{1}=[24 27 30 32]
% sts(1).revruns{2}=[46 52 62 65 70]
sts(4).sfact=[3000 2000];
sts(4).bname='bk59w37'
sts(4).initshiftind=[35 18 29]
sts(4).skipasymp=0;

sts(5).path='/oriole6/r34w20/datasum'
sts(5).matfilename='datsum.mat'
sts(5).color='k'
sts(5).ntind=1;
sts(2).contrind=[] ;
sts(5).basruns=[3 4 5 14 17];
% sts(1).revruns{1}=[24 27 30 32]
% sts(1).revruns{2}=[46 52 62 65 70]
sts(5).sfact=[2000 2000];
sts(5).bname='r34w20'
sts(5).initshiftind=[6 16]
sts(5).skipasymp=0;

% %short delay
% out1_2350ma_out2_1440ma_90msdel_20ms_200Hz

% [01;34mstim918_100msdel_70ma_20msburst[0m
% [01;34mstimon718_100msdel_20msburst_90ma[0m
% [01;34mstimon718_90msdel_10msburst_70maout1_23_out2_14[0m
% [01;34mstimon718_90msdel_10msburst_90maout1_23_out2_14[0m
% 
% 
% stim719_110msdel_10ms_80ma
% [01;34mstim719_3_100msdel_80ma_10mspulse[0m
% 
% %720
% [0m[01;34mstim720_3mspulse_200Hz_80ma_105msdel[0m
% [01;34mstimon720_5mspulse_80ma_200Hz[0m
% [01;34mstimon720_5mspulse_80ma_200Hz_105msdel[0m
% 
% [01;34mstim722_5ms_85mdel_70ma_200Hz[0m
% 
% %725
% stim725_85msdel_80ma_5ms_200Hz
% 
% %initial baseline 714_files


% bs(3).path='/oriole4/bk20bk45/'
% bs(3).matfilename='pathvals1-analdata.mat'
% bs(3).color=[0.8 0.8 0.8]
% bs(3).ntind=3;
% bs(3).contrind=[4];
% bs(3).basruns=[6 11]
% bs(3).revruns{1}=[29 31]
% bs(3).sfact=[2000 1000];
% 
% bs(4).path='/oriole4/bk20bk45/'
% bs(4).matfilename='pathvals1-analdata.mat'
% bs(4).color=[0.6 0.6 0.8]
% bs(4).ntind=1;
% bs(4).contrind=[4];
% bs(4).asympvl=21.3;
% bs(4).basruns=[6 11]
% bs(4).revruns{1}=[];
% bs(4).sfact=[1000 1000];
% 
% bs(5).path='/oriole6/bk61w42/'
% bs(5).matfilename='pathvals1-analdata.mat'
% bs(5).color=[0.6 0.6 0.8]
% bs(5).ntind=1;
% bs(5).contrind=[3];
% bs(5).asympvl=3.96;
% bs(5).basruns=[8 11 14]
% bs(5).sfact=[3000 1000];
% 
% bs(6).path='/oriole6/bk28w6/'
% bs(6).matfilename='pathvals1-analdata.mat'
% bs(6).color=[0.6 0.6 0.8]
% bs(6).ntind=1;
% bs(6).contrind=[2];
% bs(6).basruns=[24];
% bs(6).sfact=[2000 2000];
% 
% bs(7).path='/oriole/pk32bk28/'
% bs(7).matfilename='pathvals3-analdata.mat'
% bs(7).color=[0.6 0.6 0.8]
% bs(7).ntind=1;
% bs(7).contrind=[2];
% bs(7).shiftride{1}=[20];
% bs(7).basruns=[3]
% bs(7).sfact=[3000 1000];
% 
% bs(8).path='/oriole6/bk61w42/'
% bs(8).matfilename='pathvals3-analdata.mat'
% bs(8).color=[0.6 0.6 0.8]
% bs(8).ntind=1;
% bs(8).contrind=[2];
% % bs(8).asympvl=4.3
% bs(8).basruns=[2]
% bs(8).sfact=[3000 1000];
% 
% 
% 
% bs(9).path='/oriole/bk63w43/'
% bs(9).matfilename='pathvals1-analdata.mat'
% bs(9).color=[0.6 0.6 0.8]
% bs(9).ntind=1;
% bs(9).contrind=[3];
% bs(9).basruns=[5 8 11]
% bs(9).sfact=3000;
% 
% 
% 
% bs(10).path='/oriole5/pu57w52/'
% bs(10).matfilename='pathvals1-analdata.mat'
% bs(10).color=[0.6 0.6 0.8]
% bs(10).ntind=1;
% bs(10).contrind=[2];
% bs(10).asympvl=5.8;
% bs(10).basruns=[9 11]
% bs(10).sfact=3000;
% 
% bs(11).path='/oriole5/pu56w26/'
% bs(11).matfilename='pathvals1-analdata.mat'
% bs(11).color=[0.6 0.6 0.8]
% bs(11).ntind=1;
% bs(11).contrind=[2];
% bs(11).basruns=[2]
% bs(11).sfact=3000;
% 
% 
% 
% % 
% % 
% % bs(10).path='/oriole5/bk15bk14/'
% % bs().matfilename='pathvals3-analdata.mat'
% % bs(8).color=[0.6 0.6 0.8]
% % bs(8).ntind=2;
% % bs(8).contrind=[4];
% % bs(8).basruns=[2 5 7];
% % 
% % 
% % 
% % bs(11).path='/oriole5/bk15bk14/'
% % bs(11).matfilename='pathvals3-analdata.mat'
% % bs(11).color=[0.6 0.6 0.8]
% % bs(11).ntind=1;
% % bs(11).contrind=[4];
% % bs(11).basruns=[2 5 7]
% % 
