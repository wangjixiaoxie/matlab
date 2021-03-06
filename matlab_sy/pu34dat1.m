%designed to call pitchsyntaxanal, and then inactivanal.m
clear date
clear avls
clear ac tbinshft PRENT PSTNT muanal
clear graphvals edges muanal revanal wn
clear tmn;clear tmf;
sumpath='/oriole4/pu34/'
strcmd=['cd ' sumpath]
eval(strcmd);
if (~exist('datasum','dir'))
    !mkdir datasum
end
matfilename='pathvals1'
clear pathvl vals wn;

%this is calculated with two days of canin data.
avls.initmean{1}=6870
avls.maxmeans=[6614 7139]
avls.initsd{1}=80
avls.exsongpath=1;
avls.songind=30;
avls.songtm=[2.05 3.05]

avls.initmean{2}=4400;
avls.initsd{2}=101.8;

avls.initmean{3}=4493;
avls.initsd{3}=111.7;

avls.changenote=0;
%%%TKTK
edges{1}=[6300:40:7500]
edges{2}=[1800:60:2500]
edges{3}=[1800:60:2500]

wn(1).tmon{1}={'2008-7-12 7' '2008-7-16 7'}
wn(1).tmon{2}={'2008-7-27 7' '2008-7-30 7'}  
wn(1).tmon{3}={'2008-9-19 7' }
wn(1).tmlist{1}={'2008-7-12 7' '2008-7-22 7' '2008-7-22 11' '2008-7-22 17' };
wn(1).ptlist{1}=[6640 6680 6720 6820]
wn(1).tmoff{1}={'2008-7-16 7' '2008-7-22 7'}
wn(1).tmoff{2}={'2008-7-30 7' '2008-8-3 7'} 
wn(1).tmoff{3}={'2008-9-26 7'}
wn(2).tmon=wn(1).tmon;
wn(2).tmoff=wn(1).tmoff;

wn(3).tmon=wn(1).tmon;
wn(3).tmoff=wn(1).tmoff;

wnrev(1).tmon{1}={'2008-7-22 7'}
wnrev(1).tmoff{1}={'2008-7-27 7'}
wnrev(1).tmon{2}={'2008-8-3 7'}
wnrev(1).tmoff{2}={'2008-8-7 7'}
wnrev(2).tmon=wnrev(1).tmon;
wnrev(2).tmoff=wnrev(1).tmoff;
wnrev(3).tmon=wnrev(2).tmon;
wnrev(3).tmoff=wnrev(2).tmoff;
graphvals.wnrev=wnrev;


pathvl{1}='/oriole4/pu34/probein/'
catchvl{1}='batch10.keep.rand';tmn{1}=['7'];tmf{1}=['09:15'];

pathvl{2}='/oriole4/pu34/lid710/'
catchvl{2}='batchcomb';tmn{2}=['09:15'];tmf{2}=['12:30'];

pathvl{3}='/oriole4/pu34/ac710/'
catchvl{3}='batch.keep';tmn{3}=['12:30'];tmf{3}=['16:30'];

pathvl{4}='/oriole4/pu34/ac710-2/'
catchvl{4}='batch11.keep.rand';tmn{4}=['7'];tmf{4}=['09:01'];

pathvl{5}='/oriole4/pu34/lid711/'
catchvl{5}='batchcomb';tmn{5}=['09:01'];tmf{5}=['12:11'];

pathvl{6}='/oriole4/pu34/ac711/'
catchvl{6}='batch';tmn{6}=['12:11'];tmf{6}=['21:00'];

pathvl{7}='/oriole4/pu34/wnonpu34/'
catchvl{7}='batch14c.keep.rand';tmn{7}=['7'];tmf{7}=['15:30'];

pathvl{8}='/oriole4/pu34/lid714/'
catchvl{8}='batchcomb';tmn{8}=['15:30'];tmf{8}=['18:30'];

% pathvl{9}='/oriole4/pu34/ac714/'
% catchvl{9}='batch.catch';tmn{9}=['7'];tmf{9}=['21'];

pathvl{9}='/oriole4/pu34/ac714/'
catchvl{9}='batch16comb';tmn{9}=['7'];tmf{9}=['13:08'];

pathvl{10}='/oriole4/pu34/lid716/'
catchvl{10}='batchcomb';tmn{10}=['13:08'];tmf{10}=['15:44'];

pathvl{11}='/oriole4/pu34/ac714/'
catchvl{11}='batchcomb';tmn{11}=['18:30'];tmf{11}=['21'];

pathvl{12}='/oriole4/pu34/ac716/'
catchvl{12}='batch16.catch.keep';tmn{12}=['15:44'];tmf{12}=['21'];

pathvl{13}='/oriole4/pu34/ac716/'
catchvl{13}='batch17.catch.keep';tmn{13}=['7'];tmf{13}=['11:00'];


pathvl{14}='/oriole4/pu34/lid717/'
catchvl{14}='batch.catch.keep';tmn{14}=['11:00'];tmf{14}=['14:32'];

pathvl{15}='/oriole4/pu34/ac717/'
catchvl{15}='batch17.catch.keep';tmn{15}=['14:32'];tmf{15}=['21'];


pathvl{16}='/oriole4/pu34/ac717/'
catchvl{16}='batch19.catch.keep';tmn{16}=['7'];tmf{16}=['11:08'];


pathvl{17}='/oriole4/pu34/lid719/'
catchvl{17}='batch.catch.keep';tmn{17}=['11:08'];tmf{17}=['14:44'];

pathvl{18}='/oriole4/pu34/canin/'
catchvl{18}='batch06.keep.rand';tmn{18}=['7'];tmf{18}=['21'];

pathvl{19}='/oriole4/pu34/canin/'
catchvl{19}='batch07.keep.rand';tmn{19}=['7'];tmf{19}=['21'];

pathvl{20}='/oriole4/pu34/ac719/'
catchvl{20}='batch21.keep.catch';tmn{20}=['7'];tmf{20}=['10:17'];

pathvl{21}='/oriole4/pu34/lid721/'
catchvl{21}='batch.catch.keep';tmn{21}=['10:17'];tmf{21}=['13:58'];

pathvl{22}='/oriole4/pu34/ac721/'
catchvl{22}='batch21.catch';tmn{22}=['13:58'];tmf{22}=['21'];

pathvl{23}='/oriole4/pu34/wnswitch721/'
catchvl{23}='batch22.catch.keep';tmn{23}=['7'];tmf{23}=['13:35'];


pathvl{24}='/oriole4/pu34/lid722/'
catchvl{24}='batch22.keep.catch';tmn{24}=['13:35'];tmf{24}=['16:47'];

pathvl{25}='/oriole4/pu34/ac722/'
catchvl{25}='batch22.catch.keep';tmn{25}=['16:47'];tmf{25}=['21'];

pathvl{26}='/oriole4/pu34/ac722/'
catchvl{26}='batch23.catch.keep';tmn{26}=['7'];tmf{26}=['12:46'];

pathvl{27}='/oriole4/pu34/lid723/'
catchvl{27}='batch.catch.keep';tmn{27}=['12:46'];tmf{27}=['16:11'];

pathvl{28}='/oriole4/pu34/ac723/'
catchvl{28}='batch23.catch.keep';tmn{28}=['16:11'];tmf{28}=['21'];

pathvl{29}='/oriole4/pu34/ac723/'
catchvl{29}='batch25.keep.rand';tmn{29}=['7'];tmf{29}=['10:39'];

pathvl{30}='/oriole4/pu34/lid725/'
catchvl{30}='batch.keep';tmn{30}=['10:39'];tmf{30}=['14:45'];

pathvl{31}='/oriole4/pu34/ac725/'
catchvl{31}='batch25.keep.rand';tmn{31}=['14:45'];tmf{31}=['21'];

pathvl{32}='/oriole4/pu34/lid726/'
catchvl{32}='batchcomb';tmn{32}=['10:12'];tmf{32}=['13:21'];

pathvl{33}='/oriole4/pu34/ac725/'
catchvl{33}='batch26.keep.rand';tmn{33}=['7'];tmf{33}=['10:12'];

pathvl{34}='/oriole4/pu34/ac726/'
catchvl{34}='batch.keep.rand';tmn{34}=['13:21'];tmf{34}=['21'];

pathvl{35}='/oriole4/pu34/wnon2UP/'
catchvl{35}='batch28.catch.keep';tmn{35}=['7'];tmf{35}=['09:55'];

pathvl{36}='/oriole4/pu34/lid728/'
catchvl{36}='batch.catch.keep';tmn{36}=['09:55'];tmf{36}=['13:47'];

pathvl{37}='/oriole4/pu34/ac728/'
catchvl{37}='batch.catch.keep';tmn{37}=['13:47'];tmf{37}=['21:00'];

pathvl{38}='/oriole4/pu34/ac728/'
catchvl{38}='batch30.catch.keep';tmn{38}=['7'];tmf{38}=['12:04'];

pathvl{39}='/oriole4/pu34/lido730/'
catchvl{39}='batch.catch.keep';tmn{39}=['12:04'];tmf{39}=['16:17'];
 
pathvl{40}='/oriole4/pu34/ac730/'
catchvl{40}='batch30.keep';tmn{40}=['16:17'];tmf{40}=['21'];

pathvl{41}='/oriole4/pu34/ac730/'
catchvl{41}='batch01.catch.keep';tmn{41}=['7'];tmf{41}=['9:47'];

pathvl{42}='/oriole4/pu34/lid801/'
catchvl{42}='batch.keep.catch';tmn{42}=['9:47'];tmf{42}=['13:39'];

pathvl{43}='/oriole4/pu34/lid802/'
catchvl{43}='batch.keep.catch';tmn{43}=['9:29'];tmf{43}=['13:27'];

pathvl{44}='/oriole4/pu34/ac801/'
catchvl{44}='batch02.catch.keep';tmn{44}=['7'];tmf{44}=['09:30'];

pathvl{45}='/oriole4/pu34/ac802/'
catchvl{45}='batchcomb';tmn{45}=['13:27'];tmf{45}=['21'];

pathvl{46}='/oriole4/pu34/wnswitch2lidon/'
catchvl{46}='batch03lid.catch.keep';tmn{46}=['11:55'];tmf{46}=['16:08'];

pathvl{47}='/oriole4/pu34/wnswitch2lidon/'
catchvl{47}='batch04ac1.catch.keep';tmn{47}=['7'];tmf{47}=['11:30'];

pathvl{48}='/oriole4/pu34/wnswitch2lidon/'
catchvl{48}='batch04lid';tmn{48}=['11:30'];tmf{48}=['15:10'];

pathvl{49}='/oriole4/pu34/wnswitch2lidon/'
catchvl{49}='batch04ac2.catch.keep';tmn{49}=['15:10'];tmf{49}=['21'];

pathvl{50}='/oriole4/pu34/wnswitch2lidon/'
catchvl{50}='batchcomb';tmn{50}=['16:08'];tmf{50}=['21'];

pathvl{51}='/oriole4/pu34/wnswitch2lidon/'
catchvl{51}='batch05ac1.catch.keep';tmn{51}=['7'];tmf{51}=['11:47'];

pathvl{52}='/oriole4/pu34/wnswitch2lidon/'
catchvl{52}='batch05lid';tmn{52}=['11:47'];tmf{52}=['15:18'];

pathvl{53}='/oriole4/pu34/wnswitch2lidon/'
catchvl{53}='batch05ac2.catch.keep';tmn{53}=['15:18'];tmf{53}=['21'];

pathvl{54}='/oriole4/pu34/wnswitch2lidon/'
catchvl{54}='batch06.catch.keep';tmn{54}=['7'];tmf{54}=['21'];




pathvl{55}='/oriole4/pu34/wnswitch2/'
catchvl{55}='batch03.catch.keeplim';tmn{55}=['7'];tmf{55}=['11:55'];

pathvl{56}='/oriole4/pu34/ac801/'
catchvl{56}='batch01.catch.keep';tmn{56}=['13:39'];tmf{56}=['21:00'];

pathvl{57}='/oriole4/pu34/ac719/'
catchvl{57}='batch19.catch.keep';tmn{57}=['14:44'];tmf{57}=['21:00'];

pathvl{58}='/oriole4/pu34/probin3/'
catchvl{58}='batch.keep.rand';tmn{1}=['58'];tmf{58}=['8:38'];

pathvl{59}='/oriole4/pu34/lid913/'
catchvl{59}='batch.keep';tmn{59}=['8:38'];tmf{59}=['12:13'];

pathvl{60}='/oriole4/pu34/ac913/'
catchvl{60}='batch13.keep.rand';tmn{60}=['12:13'];tmf{60}=['21'];

pathvl{61}='/oriole4/pu34/ac913/'
catchvl{61}='batch14.keep.rand';tmn{61}=['7'];tmf{61}=['13:45'];

pathvl{62}='/oriole4/pu34/lid914/'
catchvl{62}='batch.keep';tmn{62}=['13:45'];tmf{62}=['17:15'];

pathvl{63}='/oriole4/pu34/ac914/'
catchvl{63}='batch14.keep';tmn{63}=['17:15'];tmf{63}=['21'];

pathvl{64}='/oriole4/pu34/ac914/'
catchvl{64}='batch15.keep.rand';tmn{64}=['7'];tmf{64}=['8:57'];

pathvl{65}='/oriole4/pu34/lid915/'
catchvl{65}='batch.keep';tmn{65}=['8:57'];tmf{65}=['12:19'];

pathvl{66}='/oriole4/pu34/ac915/'
catchvl{66}='batch.keep';tmn{66}=['12:19'];tmf{66}=['15:32'];

pathvl{67}='/oriole4/pu34/500mu915/'
catchvl{67}='batch.keep';tmn{67}=['15:32'];tmf{67}=['18:19'];

pathvl{68}='/oriole4/pu34/ac915-2/'
catchvl{68}='batch15.keep';tmn{68}=['18:19'];tmf{68}=['21'];

pathvl{69}='/oriole4/pu34/ac915-2/'
catchvl{69}='batch16.keep.rand';tmn{69}=['7'];tmf{69}=['8:23'];

pathvl{70}='/oriole4/pu34/500mu916/'
catchvl{70}='batch.keep.rand';tmn{70}=['8:23'];tmf{70}=['11:50'];

pathvl{71}='/oriole4/pu34/ac916/'
catchvl{71}='batch.keep.rand';tmn{71}=['11:50'];tmf{71}=['21'];

pathvl{72}='/oriole4/pu34/wnon916/'
catchvl{72}='batch19';tmn{72}=['7'];tmf{72}=['21'];

pathvl{73}='/oriole4/pu34/ac920/'
catchvl{73}='batch22.keep.catch';tmn{73}=['7'];tmf{73}=['11:29'];

pathvl{74}='/oriole4/pu34/mu922/'
catchvl{74}='batch22.keep.catch';tmn{74}=['11:29'];tmf{74}=['14:59'];

pathvl{75}='/oriole4/pu34/ac922/'
catchvl{75}='batch22.keep.catch';tmn{75}=['14:59'];tmf{75}=['21'];

pathvl{76}='/oriole4/pu34/wnon916/'
catchvl{76}='batch20.catch.keep';tmn{76}=['7'];tmf{76}=['11:51'];

pathvl{77}='/oriole4/pu34/500mu920/'
catchvl{77}='batch.catch.keep';tmn{77}=['11:51'];tmf{77}=['15:26'];

pathvl{78}='/oriole4/pu34/ac920/'
catchvl{78}='batchcomb';tmn{78}=['15:26'];tmf{78}=['21'];

pathvl{79}='/oriole4/pu34/wnonpu34/'
catchvl{79}='batch12.catch.keep';tmn{79}=['7'];tmf{79}=['21'];

pathvl{80}='/oriole4/pu34/wnonpu34/'
catchvl{80}='batch13.catch.keep';tmn{80}=['7'];tmf{80}=['21'];

pathvl{81}='/oriole4/pu34/ac714/'
catchvl{81}='batch15.catch.keep';tmn{81}=['7'];tmf{81}=['21'];

pathvl{82}='/oriole4/pu34/wnon2UP/'
catchvl{82}='batch27.catch.keep';tmn{82}=['7'];tmf{82}=['21'];

pathvl{83}='/oriole4/pu34/ac728/'
catchvl{83}='batch29.catch.keep';tmn{83}=['7'];tmf{83}=['21'];

%adding extra data 1.12.2010

pathvl{84}='/oriole4/pu34/ac922/'
catchvl{84}='batch24.keep.catch';tmn{84}=['7'];tmf{84}=['9:57'];

pathvl{85}='/oriole4/pu34/500mu924/'
catchvl{85}='batchcomb';tmn{85}=['9:57'];tmf{85}=['13:51'];

pathvl{86}='/oriole4/pu34/ac924/'
catchvl{86}='batchcomb';tmn{86}=['13:51'];tmf{86}=['21'];


pathvl{87}='/oriole4/pu34/ac924/'
catchvl{87}='batch25.keep.rand';tmn{87}=['7'];tmf{87}=['8:37'];

pathvl{88}='/oriole4/pu34/500mu925/'
catchvl{88}='batch.keep';tmn{88}=['8:37'];tmf{88}=['12:42'];


pathvl{89}='/oriole4/pu34/iboacid925/'
catchvl{89}='batch';tmn{89}=['12:42'];tmf{89}=['21'];

pathvl{90}='/oriole4/pu34/ac925-2/'
catchvl{90}='batch26.keep';tmn{90}=['7'];tmf{90}=['10'];

pathvl{91}='/oriole4/pu34/ac926/'
catchvl{91}='batch27.catch.keep';tmn{91}=['7'];tmf{91}=['21'];


pathvl{92}='/oriole4/pu34/ac926/'
catchvl{92}='batch28.catch.keep';tmn{92}=['7'];tmf{92}=['21'];


pathvl{93}='/oriole4/pu34/ac926/'
catchvl{93}='batch29.catch.keep';tmn{93}=['7'];tmf{93}=['21'];


pathvl{94}='/oriole4/pu34/ac926/'
catchvl{94}='batch30.catch.keep';tmn{94}=['7'];tmf{94}=['21'];


pathvl{95}='/oriole4/pu34/wnon1027/'
catchvl{95}='batch.catch.keep2';tmn{95}=['7'];tmf{95}=['21'];
% 
pathvl{96}='/oriole4/pu34/ac926/'
catchvl{96}='batch02.catch.keep.rand';tmn{96}=['7'];tmf{96}=['21'];

pathvl{97}='/oriole4/pu34/ac717/'
catchvl{97}='batch18.catch';tmn{97}=['7'];tmf{97}=['21'];


pathvl{98}='/oriole4/pu34/ac719/'
catchvl{98}='batch20.catch';tmn{98}=['7'];tmf{98}=['21'];


pathvl{99}='/oriole4/pu34/ac920/'
catchvl{99}='batch21';tmn{99}=['7'];tmf{99}=['21'];

pathvl{100}='/oriole4/pu34/ac922/'
catchvl{100}='batch23';tmn{100}=['7'];tmf{100}=['21'];


usex(1:49)=0;
ac(49)=[0];
avls.analind=[1:49]

%how to deal with dirfiles -- put all in one directory.
dirpathvl{1}='/oriole4/pk32bk28/dirfiles/'
dircatchvl{1}='batch'

dirpathvl{2}='/oriole4/pk32bk28/wnon428/'
dircatchvl{2}='batchcomb'


numnotes=2
notes='ab'

fbins={};
tbinshft{1}=0.076;
tbinshft{2}=0.008;  
tbinshft{3}=0.005;
avls.pitchalign{1}=[0]
avls.pitchalign{2}=0;
avls.pitchalign{3}=0;



NFFT(1)=512
NFFT(2)=256
NFFT(3)=512; 

fbins{1}=[6400,7400];
fbins{2}=[4000,6000];
fbins{3}=[4000,6000];

clear NT     
NT{1}='a'
NT{2}='b' 

 NT{3}='b'
avls.exclude=[0 1 1]
PRENT{1}='';PSTNT{1}='';
PRENT{2}='a';PSTNT{2}='b';
PRENT{3}='b';PSTNT{3}='b';
  %%%plotvals
  
acon=[ 1 3 4 6 7 9 11 12 13 15  16 18 19 20 22 23 25 26 28 29 31 33 34 35 37 38 40 41 44 45 47 49 50 51 53 54 55 56 57 58 60 61 63 64 68 69 71 72 73 75 76 78:83 84 86:87 89:100  ] ;  
muon=[2 5 8 10 14 17 21 24 27 30 32 36 39 42 43 46 48 52 59 62 65 66 67 70 74 77 85   ]
revanal{1}{1}=[24 27 30];
revanal{1}{2}=[46 48 52]
revanal{2}{1}=[24 27 30]
revanal{2}{2}=[46 48 52]
muanal{1}{1}=[2 5 8 10 14 17 21 ]
muanal{1}{2}=[32 36 39 42 43] 
muanal{2}{1}=[2 5 8 10 14 17 21 ]
muanal{2}{2}=[32 36 39 42 43] 


diron=[];
diron=[]
extraon=[]

clear colvals
colvals{1}=muon;
colvals{2}=acon
colvals{3}=diron;

graphvals.numcol=3
graphvals.col='rkc'
% graphvals.plttext=1
% graphvals.txtht=[3260 3200]
% graphvals.tickht=8000
% graphvals.edges{1}=[6700:80:7900]
% graphvals.edges{2}=[2000:40:2500]
% graphvals.edges{3}=[1300:40:1800]
% graphvals.colvals=colvals
% graphvals.pltpnts=1
avls.bname='pu34'
avls.diron=diron;
avls.usex=usex
avls.tmon=tmn;
avls.tmoff=tmf;
graphvals.tmn=tmn;
graphvals.tmf=tmf;
graphvals.acon=ac;
graphvals.date=date;
avls.edges=edges;
graphvals.wn=wn;

graphvals.colvals=colvals
graphvals.chunkdata=1

avls.muanal=muanal;
avls.revanal=revanal;
avls.pvls=pathvl
avls.datfile=matfilename;
avls.cvl=catchvl
avls.sumpath=sumpath
avls.mtflnm=matfilename  
avls.supanal=0
avls.NT=NT
avls.acon=acon;
avls.muon=muon;
avls.NFFT=NFFT
avls.fbins=fbins
avls.tshft=tbinshft
avls.numnotes=numnotes
avls.mkfv=[]
avls.PRENT=PRENT
avls.PSTNT=PSTNT
avls.repeatanal=[0 0 0]
avls.deadtm=.0833/24;
avls.muoffset=1.5/24;
avls.acoffset=1.5/24;
avls.pretm=4/24;


strcmd=['cd ' sumpath 'datasum']
eval(strcmd)

strcmd=['save ' matfilename '.mat avls graphvals'];
eval(strcmd);


