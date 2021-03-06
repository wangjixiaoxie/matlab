%this script is written to generate pitch contours for bk57w35
%writes data to a .mat file in directory where data is.

%YOU NEED TO SET MUINDS BELOW
%calls jc_pitchcontourFV, in ~/matlab
%      maketimevec2
clear fv fbins bt smoothamp baspath conbins avls


avls.baspath='/oriole/bk48w74/'
avls.datdir='datasum'

wnin(1).tmon{1}={'2009-8-22 7' '2009-8-24 7'}
wnin(1).tmon{2}={'2009-9-02 7' '2009-9-04 7'}
wnin(1).tmon{3}={'2009-9-26 7'}
wnin(1).tmon{4}={'2009-10-31 7'}
wnin(1).tmoff{1}={'2009-8-24 7' '2009-8-29 7'}
wnin(1).tmoff{2}={'2009-9-04 7' '2009-9-09 7'} 
wnin(1).tmoff{3}={'2009-10-13 7'}
wnin(1).tmoff{4}={'2009-11-10 7'}


wnrevin(1).tmon{1}={'2009-8-29 7'}
wnrevin(1).tmoff{1}={'2009-8-31 7'}
wnrevin(1).tmon{2}={'2009-9-09 7'}
wnrevin(1).tmoff{2}={'2009-9-12 7'}
wnrevin(1).tmon{2}={'2009-9-09 7'}
wnrevin(1).tmoff{2}={'2009-9-12 7'}
wnrevin(1).tmon{3}={'2009-10-13 7'}
wnrevin(1).tmoff{3}={'2009-10-16 7'}
wnrevin(1).tmon{4}={'2009-11-11 7'}
wnrevin(1).tmoff{4}={'2009-11-15 7'}


avls.wnin=wnin
avls.wnrevin=wnrevin



%leadin
%temptest_aug21
%temptest

avls.pvls{1}='stim14_bilat_90ma_65msdel_200Hz_run2/'
avls.cvl{1}='batch'
avls.NT{1}='a'

avls.pvls{2}='stim14_90ma_821/'
avls.cvl{2}='batch21.keep'
avls.NT{2}='a'

avls.pvls{3}='UPSHIFT/stimaug22/'
avls.cvl{3}='batch.keep'
avls.NT{3}='a'

avls.pvls{4}='stimaug23/'
avls.cvl{4}='batch'
avls.NT{4}='a'

avls.pvls{5}='stimaug24/'
avls.cvl{5}='batch'
avls.NT{5}='a'


avls.pvls{6}='stimaug25/'
avls.cvl{6}='batch25'
avls.NT{6}='a'

avls.pvls{7}='stimaug26/'
avls.cvl{7}='batch'
avls.NT{7}='a'

avls.pvls{8}='stimaug23_65msdel_5mspulse/'
avls.cvl{8}='batch'
avls.NT{8}='a'

avls.pvls{9}='stimaug23_95msdel_5mspulse/'
avls.cvl{9}='batch23'
avls.NT{9}='a'

avls.pvls{10}='stimaug27_110msdel_10ms_CORRECT/'
avls.cvl{10}='batch'
avls.NT{10}='a'
avls.del{10}=.110

avls.pvls{11}='stimaug27_95msdel_10mspulse/'
avls.cvl{11}='batch'
avls.NT{11}='a'
avls.del{11}=.095

avls.pvls{12}='stimaug27_65msdel_10mspulse/'
avls.cvl{12}='batch'
avls.NT{12}='a'
avls.del{12}=.065;

avls.pvls{13}='stimaug28/'
avls.cvl{13}='batch'
avls.NT{13}='a'
avls.del{13}=.065;

% avls.pvls{14}='temptest829-1/'
% avls.cvl{14}='batch'
% avls.NT{14}='a'
% avls.del{14}=.065;
% 
% avls.pvls{15}='revcontinaug28/'
% avls.cvl{15}='batch'
% atvls.NT{15}='a'
% avls.del{15}=.065;

% tvls.pvls{16}='revcontinaug28/'
% tvls.cvl{16}='batch29.keep'
% tvls.NT{16}='a'
% tvls.del{16}=.065;

avls.pvls{14}='stimaug29-2/'
avls.cvl{14}='batch'
avls.NT{14}='a'
avls.del{14}=.065;

% avls.pvls{17}='revcontinaug28/'
% avls.cvl{17}='batch17'
% avls.NT{17}='a'
% avls.del{17}=.065;

avls.pvls{15}='stimaug29-3/'
avls.cvl{15}='batch'
avls.NT{15}='a'
avls.del{15}=.065;

avls.pvls{16}='95msdel_10ms_aug29/'
avls.cvl{16}='batch'
avls.NT{16}='a'
avls.del{16}=.095;

avls.pvls{17}='stimaug30/'
avls.cvl{17}='batch30'
avls.NT{17}='a'
avls.del{17}=.065

avls.pvls{18}='stimaug31/'
avls.cvl{18}='batch'
avls.NT{18}='a'
avls.del{18}=.065

avls.pvls{19}='stim901/'
avls.cvl{19}='batch'
avls.NT{19}='a'
avls.del{19}=.065

avls.pvls{20}='901_shortstim/120msdel_5mspulse'
avls.cvl{20}='batch'
avls.NT{20}='a'
avls.del{20}=.12

avls.pvls{21}='901_shortstim/100msdel_5mspulse'
avls.cvl{21}='batch'
avls.NT{21}='a'
avls.del{21}=.1

avls.pvls{22}='901_shortstim/80msdel_5mspulse'
avls.cvl{22}='batch'
avls.NT{22}='a'
avls.del{22}=.08

avls.pvls{23}='stim901-2'
avls.cvl{23}='batch'
avls.NT{23}='a'
avls.del{23}=.065

avls.pvls{24}='UPSHIFT/stim902'
avls.cvl{24}='batch'
avls.NT{24}='a'
avls.del{24}=.065

avls.pvls{25}='UPSHIFT/stim903-1'
avls.cvl{25}='batch.keep'
avls.NT{25}='a'
avls.del{25}=.065

avls.pvls{26}='UPSHIFT/stim903_110msdel_10ms'
avls.cvl{26}='batch03.keep'
avls.NT{26}='a'
avls.del{26}=.110

avls.pvls{27}='UPSHIFT/stim903-2'
avls.cvl{27}='batch.keep'
avls.NT{27}='a'
avls.del{27}=.065

avls.pvls{28}='UPSHIFT/stim903-2'
avls.cvl{28}='batch.keep'
avls.NT{28}='a'
avls.del{28}=.065

avls.pvls{29}='UPSHIFT/stim904_65msdel'
avls.cvl{29}='batch.keep'
avls.NT{29}='a'
avls.del{29}=.065

avls.pvls{30}='UPSHIFT/stim904_110msdel'
avls.cvl{30}='batch.keep'
avls.NT{30}='a'
avls.del{30}=.110

avls.pvls{31}='UPSHIFT/stim907standard'
avls.cvl{31}='batch.keep'
avls.NT{31}='a'
avls.del{31}=.065

avls.pvls{32}='UPSHIFT/stim907_110msdel'
avls.cvl{32}='batch.keep'
avls.NT{32}='a'
avls.del{32}=.110

avls.pvls{33}='UPSHIFT/stim908standard'
avls.cvl{33}='batch.keep'
avls.NT{33}='a'
avls.del{33}=.065


avls.pvls{34}='REVSHIFT909/wnonrev909'
avls.cvl{34}='batch10.keep'
avls.NT{34}='a'
avls.del{34}=.065

avls.pvls{35}='REVSHIFT909/stim909-1standard'
avls.cvl{35}='batch'
avls.NT{35}='a'
avls.del{35}=.065

avls.pvls{36}='REVSHIFT909/stim909-2standard'
avls.cvl{36}='batch.keep'
avls.NT{36}='a'
avls.del{36}=.065


avls.pvls{37}='REVSHIFT909/stim909_125msdel'
avls.cvl{37}='batch.keep'
avls.NT{37}='a'
avls.del{37}=.125


avls.pvls{38}='REVSHIFT909/stim909_125-2'
avls.cvl{38}='batch.keep'
avls.NT{38}='a'
avls.del{38}=.125

avls.pvls{39}='REVSHIFT909/stim910-standard'
avls.cvl{39}='batch.keep'
avls.NT{39}='a'
avls.del{39}=.065

avls.pvls{40}='REVSHIFT909/stim910_110msdel'
avls.cvl{40}='batch.keep'
avls.NT{40}='a'
avls.del{40}=.110


avls.pvls{41}='REVSHIFT909/stim911-standard'
avls.cvl{41}='batch11.keep'
avls.NT{41}='a'
avls.del{41}=.065

avls.pvls{42}='REVSHIFT909/stim912-standard'
avls.cvl{42}='batch.keep'
avls.NT{42}='a'
avls.del{42}=.065

avls.pvls{43}='REVSHIFT909/stim913-standard'
avls.cvl{43}='batch.keep'
avls.NT{43}='a'
avls.del{43}=.065

avls.pvls{44}='stimaug27_85msdel_10mspulse/'
avls.cvl{44}='batch'
avls.NT{44}='a'
avls.del{44}=.0865;

avls.pvls{45}='DOWNSHIFT2/bas_stim_90msdel_60ms_925/'
avls.cvl{45}='batch.keep'
avls.NT{45}='a'
avls.del{45}=.090;

avls.pvls{46}='DOWNSHIFT2/bas_stim_90msdel_60ms_50MA_925/'
avls.cvl{46}='batch.keep'
avls.NT{46}='a'
avls.del{46}=.090;

avls.pvls{47}='DOWNSHIFT2/bas_stim_110msdel_10ms_925/'
avls.cvl{47}='batch.keep'
avls.NT{47}='a'
avls.del{47}=.110;

avls.pvls{48}='DOWNSHIFT2/stim928_10ms_100ma_70msdel-2/'
avls.cvl{48}='batch.keep'
avls.NT{48}='a'
avls.del{48}=.07;

avls.pvls{49}='DOWNSHIFT2/stim928_10ms_100ma_90msdel-2/'
avls.cvl{49}='batch.keep'
avls.NT{49}='a'
avls.del{49}=.09;

avls.pvls{50}='DOWNSHIFT2/stim928_10ms_100ma_110msdel-2/'
avls.cvl{50}='batch.keep'
avls.NT{50}='a'
avls.del{50}=.11;

avls.pvls{51}='DOWNSHIFT2/stim928_10ms_100ma_130msdel-2/'
avls.cvl{51}='batch.keep'
avls.NT{51}='a'
avls.del{51}=.13;

avls.pvls{52}='DOWNSHIFT2/COR_STIM929_75msdel_60ms_90ma/'
avls.cvl{52}='batch.keep'
avls.NT{52}='a'
avls.del{52}=.075;

avls.pvls{53}='stim929_100ma_90msdel_10ms/'
avls.cvl{53}='batch.keep'
avls.NT{53}='a'
avls.del{53}=.09;

avls.pvls{54}='DOWNSHIFT2/stim1001_75msdel_60ms_90ma/'
avls.cvl{54}='batch.keep'
avls.NT{54}='a'
avls.del{54}=.075;


avls.pvls{55}='DOWNSHIFT2/stim926-3-75msdel/'
avls.cvl{55}='batch.keep'
avls.NT{55}='a'
avls.del{55}=.075;

%broken wire
avls.pvls{56}='DOWNSHIFT2/stim1005_90ma_75msdel_60ms/'
avls.cvl{56}='batch.keep'
avls.NT{56}='a'
avls.del{56}=.075;

avls.pvls{57}='DOWNSHIFT2/stim1006_COR2_60ms/'
avls.cvl{57}='batch.keep'
avls.NT{57}='a'
avls.del{57}=.075;

avls.pvls{58}='DOWNSHIFT2/stim1008_60ms_75msdel_90ma/'
avls.cvl{58}='batch.keep'
avls.NT{58}='a'
avls.del{58}=.090;


avls.pvls{59}='DOWNSHIFT2/stim1012_60ms_75msdel/'
avls.cvl{59}='batch.keep'
avls.NT{59}='a'
avls.del{59}=.075;


avls.pvls{60}='DOWNSHIFT2/stim1009_90msdel_10ms/'
avls.cvl{60}='batch.keep'
avls.NT{60}='a'
avls.del{60}=.090;

avls.pvls{61}='REVSHIFT1012/60ms_75msdel_1012/'
avls.cvl{61}='batch.keep'
avls.NT{61}='a'
avls.del{61}=.090;

avls.pvls{62}='REVSHIFT1012/60ms_75msdel_1013_2/'
avls.cvl{62}='batch.keep.rand'
avls.NT{62}='a'
avls.del{62}=.090;

avls.pvls{63}='REVSHIFT1012/75msdel_60ms_1015/'
avls.cvl{63}='batch.keep'
avls.NT{63}='a'
avls.del{63}=.090;

avls.pvls{64}='stim1028/'
avls.cvl{64}='batch.keep'
avls.NT{64}='a'
avls.del{64}=.090;

avls.pvls{65}='stim1030_10ms_90msdel/'
avls.cvl{65}='batch.keep'
avls.NT{65}='a'
avls.del{65}=.090;

avls.pvls{66}='BASELINE4/stim1030_out123_out214/'
avls.cvl{66}='batch.keep'
avls.NT{66}='a'
avls.del{66}=.090;

avls.pvls{67}='BASELINE4/stim_60ms_75msdel_out123_out214/'
avls.cvl{67}='batch.keep'
avls.NT{67}='a'
avls.del{67}=.090;

avls.pvls{68}='BASELINE4/stim14_90ma_60ms_75msdelC/'
avls.cvl{68}='batch.keep'
avls.NT{68}='a'
avls.del{68}=.090;

avls.pvls{69}='BASELINE4/stim14_75msdel_60ms/'
avls.cvl{69}='batch.keep'
avls.NT{69}='a'
avls.del{69}=.075;

avls.pvls{70}='BASELINE4/stim14_90ma_1101-b_60ms/'
avls.cvl{70}='batch.keep'
avls.NT{70}='a'
avls.del{70}=.075;

avls.pvls{71}='BASELINE4/SHORTPULSE/stim110del_10msCOR/'
avls.cvl{71}='batch.keep'
avls.NT{71}='a'
avls.del{71}=.110

avls.pvls{72}='BASELINE4/SHORTPULSE/stim70del_10msCOR/'
avls.cvl{72}='batch.keep'
avls.NT{72}='a'
avls.del{72}=.070;


avls.pvls{73}='BASELINE4/SHORTPULSE/stim90del_10ms/'
avls.cvl{73}='batch.keep'
avls.NT{73}='a'
avls.del{73}=.090;

avls.pvls{74}='BASELINE4/SHORTPULSE/stim130del_10ms/'
avls.cvl{74}='batch.keep'
avls.NT{74}='a'
avls.del{74}=.130;

avls.pvls{75}='BASELINE4/stim1102_75msdel_60ms/'
avls.cvl{75}='batch.keep'
avls.NT{75}='a'
avls.del{75}=.075;

avls.pvls{76}='BASELINE4/stim90del_10ms_1103/'
avls.cvl{76}='batch.keep'
avls.NT{76}='a'
avls.del{76}=.090;

avls.pvls{77}='BASELINE4/stim14_60ms_75msdel_1104/'
avls.cvl{77}='batch.keep'
avls.NT{77}='a'
avls.del{77}=.090;

avls.pvls{78}='BASELINE4/stim1106_60ms_75msdel/'
avls.cvl{78}='batch.keep'
avls.NT{78}='a'
avls.del{78}=.075;

avls.pvls{79}='BASELINE4/stim1107_60ms_75msdel/'
avls.cvl{79}='batch.keep'
avls.NT{79}='a'
avls.del{79}=.075;

avls.pvls{80}='BASELINE4/stim1108_60ms_75msdel/'
avls.cvl{80}='batch.keep'
avls.NT{80}='a'
avls.del{80}=.075;

avls.pvls{81}='BASELINE4/stim1109_60ms_75msdel/'
avls.cvl{81}='batch.keep'
avls.NT{81}='a'
avls.del{81}=.075;



avls.pvls{82}='BASELINE4/stim1111REAL_60ms_75msdel/'
avls.cvl{82}='batch.keep'
avls.NT{82}='a'
avls.del{82}=.075;


avls.pvls{83}='BASELINE4/stim1111_2REAL/'
avls.cvl{83}='batch'
avls.NT{83}='a'
avls.del{83}=.075;


avls.pvls{84}='BASELINE4/stim1112/'
avls.cvl{84}='batch12.keep'
avls.NT{84}='a'
avls.del{84}=.075;

avls.pvls{85}='BASELINE4/stim1114/'
avls.cvl{85}='batch14.keep'
avls.NT{85}='a'
avls.del{85}=.075;

avls.pvls{86}='BASELINE4/stim1115-b/'
avls.cvl{86}='batch15.keep.rand'
avls.NT{86}='a'
avls.del{86}=.075;

avls.pvls{87}='BASELINE4/stim1116/'
avls.cvl{87}='batch.keep.rand'
avls.NT{87}='a'
avls.del{87}=.075;


for ii=1:8
    
    avls.del{ii}=.065
end
for ii=9
    avls.del{ii}=.095;
end
%do pitch analysis on these batch inds.

avls.mkfv=[]

    avls.con_tbinshft=.05;
    %shift for pitch analysis
    avls.pt_tbinshft=.077;
    %pitch
    avls.pt_NFFT=512
    %contours
    avls.con_NFFT=4096;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
  %shift for contours
   
    avls.fbins=[6100 8000]
%     avls.conbins=[2400 2800];
    
    
    avls.basruns=[]
    avls.contanal=1
avls.contfv=[]
avls.ptfv=[]
avls.analfv=[1:87]
avls.catchstimfv=[1:87]
avls.contnum=3;
    avls.con_tbinshft=.05;
    %shift for pitch analysis
    avls.pt_tbinshft=.079;
    %pitch
    avls.pt_NFFT=512
    %contours
    avls.con_NFFT=4096;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
  %shift for contours
    avls.mnbas=[7172]
    avls.stdbas=[118]
    avls.basruns=1;
    
    avls.SOUNDSTR='wav'
    avls.STIMSTR='rig'
    avls.STIMBNDS=[200 0]
     avls.HST_EDGES=[6300:50:8000]
%     avls.STAN_RUNS=[1:7 13:15 17:19 23:25 27:36 39 41:52 54:56]
    avls.STAN_RUNS=[1:7 13:15 17:19 23:25 28:33 35:36 39 41:43 45 52 54:55 57:63 65:68 70 75 77:87]
    avls.REMOVEOUTLIERS=1;
    avls.NUMSTD=2;
    conbinsmid=[30]	
    conbinshi=[73:87]
    conbinsnorm=[1:29 31:72]
for ii=1:length(conbinshi)
    avls.conbins{conbinshi(ii)}= [2400 2800]
end

for ii=1:length(conbinsmid)
    avls.conbins{conbinsmid(ii)}= [2270 2800]
end

for ii=1:length(conbinsnorm)
    avls.conbins{conbinsnorm(ii)}=[2100 2800]
end

    
    
    
    
    