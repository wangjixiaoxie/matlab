%rewritten 1.22.2010 for purposes of adding pre and post on reverse run.
%this script is written to generate pitch contours for bk61w42.
clear fv fbins bt smoothamp baspath conbins avls


avls.baspath='/oriole6/bk59w37redo/'
avls.datdir='datasum'

clear wnin
clear wnrevin

wnin(1).tmon{1}={'2009-7-15 7' '2009-7-16 7'}
wnin(1).tmon{2}={'2009-7-23 7' '2009-7-25 7'}  
wnin(1).tmoff{1}={'2009-7-16 7' '2009-7-19 7'}
wnin(1).tmoff{2}={'2009-7-25 7' '2009-7-29 7'} 

wnrevin(1).tmon{1}={'2009-7-19 7'}
wnrevin(1).tmoff{1}={'2009-7-22 7'}
wnrevin(1).tmon{2}={'2009-7-29 7'}
wnrevin(1).tmoff{2}={'2009-8-03 21'}

avls.wnin=wnin
avls.wnrevin=wnrevin
% 
% 
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

avls.pvls{1}='714/out1_23_40ma_out2_1450ma_200Hz_70msdel/'
avls.cvl{1}='batch'
avls.NT{1}='a'

avls.pvls{2}='714/out1_23_40ma_out2_1450ma_200Hz_70msdel/'
avls.cvl{2}='batch'
avls.NT{2}='a'

avls.pvls{3}='714/out1_23_40ma_out2_1450ma_200Hz_70msdel/'
avls.cvl{3}='batch'
avls.NT{3}='a'

%716files -- first day of pitchshift
avls.pvls{4}='716/stimon716_35maboth'
avls.cvl{4}='batch'
avls.NT{4}='a'

avls.pvls{5}='716/stimon716_40maout1_35maout2'
avls.cvl{5}='batch'
avls.NT{5}='a'

avls.pvls{6}='716/stimon716_40maout1_35maout2'
avls.cvl{6}='batch'
avls.NT{6}='a'

%717files  -- first day of pitchshift
avls.pvls{7}='717/stim71713out2_14out1_50ma'
avls.cvl{7}='batch'
avls.NT{7}='a'

avls.pvls{8}='717/stim71713out2_14out1_50ma'
avls.cvl{8}='batch'
avls.NT{8}='a'


% avls.pvls{8}='717/stimon717_13out2_24_out1'
% avls.cvl{8}='batch'
% avls.NT{8}='a'

avls.pvls{9}='717/stimon717_23out2_14out1_40ma_60ms_75del_200Hz'
avls.cvl{9}='batch'
avls.NT{9}='a'
%718files 

avls.pvls{10}='718/stim718_400Hz_40ms_75msdel'
avls.cvl{10}='batch'
avls.NT{10}='a'
avls.fmod{10}=[2350 2800]

avls.pvls{11}='718/stimon718_40ms_75del_out2_14_out1_23'
avls.cvl{11}='batch'
avls.NT{11}='a'
avls.del{11}=.075
avls.fmod{11}=[2350 2800]

%719 files
avls.pvls{12}='719/stim719_overnight_75msdel_60mspulse_200Hz'
avls.cvl{12}='batch1'
avls.NT{12}='a'
avls.del{12}=.075

avls.pvls{13}='719/stimprobe719'
avls.cvl{13}='batch'
avls.NT{13}='a'

avls.pvls{14}='/719/stim719-2/'
avls.cvl{14}='batch'
avls.NT{14}='a'

%720 files
avls.pvls{15}='720/stimon720_60ms_200Hz'
avls.cvl{15}='batch'
avls.NT{15}='a'

%721files
avls.pvls{16}='721/stim721_STIMONCOR_75del_200Hz_40ma_60ms'
avls.cvl{16}='batch'
avls.NT{16}='a'

%722files
avls.pvls{17}='722/stim722_40ms_40ma_75msdel'
avls.cvl{17}='batch'
avls.NT{17}='a'

avls.pvls{18}='722/stimon_standard_60ms_200Hz_75del_722'
avls.cvl{18}='batch'
avls.NT{18}='a'

%723files
avls.pvls{19}='723/stim723_75msdel_40ms_200Hz_40ma'
avls.cvl{19}='batch'
avls.NT{19}='a'

%724files???what happened 7/24
avls.pvls{20}='724/bilatstim724_45ma_40ms_75msdel'
avls.cvl{20}='batch'
avls.NT{20}='a'

%725 files
avls.pvls{21}='725/stim725_50ma_40ms_200Hz'
avls.cvl{21}='batch'
avls.NT{21}='a'
avls.fmod{21}=[2050 2600]
%726files
avls.pvls{22}='726/stim726_50ma_40ms_75del'
avls.cvl{22}='batch'
avls.NT{22}='a'

%727files
avls.pvls{23}='727/stim727_standard_50ma'
avls.cvl{23}='batch'
avls.NT{23}='a'

avls.pvls{24}='728/stim728_standard'
avls.cvl{24}='batch'
avls.NT{24}='a'

avls.pvls{25}='730/stim730'
avls.cvl{25}='batch'
avls.NT{25}='a'

avls.pvls{26}='731/stim731'
avls.cvl{26}='batch.rand'
avls.NT{26}='a'

avls.pvls{27}='801/stim801'
avls.cvl{27}='batch'
avls.NT{27}='a'

avls.pvls{28}='802/stim802'
avls.cvl{28}='batch'
avls.NT{28}='a'

avls.pvls{29}='718/71890msdelcomb'
avls.cvl{29}='batch'
avls.NT{29}='a'
avls.del{29}=.090;
avls.fmod{29}=[2350 2800]

avls.pvls{30}='718/718100msdelcomb'
avls.cvl{30}='batch'
avls.NT{30}='a'
avls.del{30}=.100
avls.fmod{30}=[2350 2800]

avls.pvls{31}='725/stim725_85msdel_80ma_5ms_200Hz'
avls.cvl{31}='batch'
avls.NT{31}='a'
avls.del{31}=.085
avls.fmod{31}=[2050 2600]


avls.pvls{32}='wnrev719'
avls.cvl{32}='batch19'
avls.NT{32}='a'
avls.del{32}=.085

avls.NOSTIM(32)=1;

avls.pvls{33}='wnrev719'
avls.cvl{33}='batch20'
avls.NT{33}='a'
avls.del{33}=.085

avls.NOSTIM(33)=1;

avls.pvls{34}='stimoff716_wnon'
avls.cvl{34}='batch18'
avls.NT{34}='a'
avls.del{34}=.085
avls.NOSTIM(34)=1;
avls.pvls{35}=avls.pvls{12}
avls.cvl{35}='batch2'
avls.NT{35}='a'
avls.del{35}=avls.del{12}

avls.NOSTIM(35)=0;


for ii=1:28
    avls.del{ii}=.075
end

avls.basruns=[1 2 3]
avls.contanal=1;
avls.contfv=[]
avls.ptfv=[]
avls.analfv=[1 4 9:35]
avls.catchstimfv=[1 4 9:35]
avls.STAN_RUNS=[1 4 9:19 21:24 27:28 32:35]
avls.contnum=3;
avls.REMOVEOUTLIERS=1;
avls.NUMSTD=2;
    avls.con_tbinshft=.05;
    %shift for pitch analysis
    avls.pt_tbinshft=.079;
    %pitch
    avls.pt_NFFT=512
    %contours
    avls.con_NFFT=4096;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
  %shift for contours
   
    avls.fbins=[6100 8000]
 
    avls.mnbas=[6981]
    avls.stdbas=[106]
    avls.basruns=1;
    
    avls.SOUNDSTR='und'
    avls.STIMSTR='im'
    avls.STIMBNDS=[200 0]
     avls.HST_EDGES=[6300:50:8000]
     
conbinshi=[10 11 21 29 30 31]
conbinsnorm=[1:9 12:20 22:28 32:35]
for ii=1:length(conbinshi)
    avls.conbins{conbinshi(ii)}=avls.fmod{conbinshi(ii)}
end
for ii=1:length(conbinsnorm)
    avls.conbins{conbinsnorm(ii)}=[2150 2800]
end
% REPSOUNDSTR='und'
% REPSTIMSTR='im'
% 
% REPINDS=[]



