%designed to call pitchsyntaxanal, and then inactivanal.m
clear date
clear avls
clear ac tbinshft PRENT PSTNT muanal
clear graphvals edges muanal revanal
clear tmn;clear tmf;
sumpath='/oriole4/pu34/'
strcmd=['cd ' sumpath]
eval(strcmd);
if (~exist('datasum','dir'))
    !mkdir datasum
end
matfilename='pathvals2'
clear pathvl vals wn;

%this is calculated with two days of canin data.
avls.initmean{1}=6874
avls.maxmeans=[6614 7139]
avls.initsd{1}=78.1
avls.exsongpath=1;
avls.songind=30;
avls.songtm=[2.05 3.05]

avls.initmean{2}=2260;
avls.initsd{2}=100;

avls.initmean{3}=2260;
avls.initsd{3}=100;

%%%TKTK
wn{1}.freqlo=[6000 7280]
wn{2}.freqhi=[7430 6000]
edges{1}=[6300:65:7500]
edges{2}=[1800:80:2500]
edges{3}=[1800:80:2500]

wn{1}.tmon={'2008-4-28 12:00' '2008-05-13 7:00'}
wn{1}.tmoff={'2008-5-13 7:00' '2008-05-16 21:00'}

pathvl{1}='/oriole4/pu34/probin3/'
catchvl{1}='batch.keep.rand';tmn{1}=['7'];tmf{1}=['8:38'];

pathvl{2}='/oriole4/pu34/lid913/'
catchvl{2}='batch.keep';tmn{2}=['8:38'];tmf{2}=['12:13'];

pathvl{3}='/oriole4/pu34/ac913/'
catchvl{3}='batch13.keep.rand';tmn{3}=['12:13'];tmf{3}=['21'];

pathvl{4}='/oriole4/pu34/ac913/'
catchvl{4}='batch14.keep.rand';tmn{4}=['7'];tmf{4}=['13:45'];

pathvl{5}='/oriole4/pu34/lid914/'
catchvl{5}='batch.keep';tmn{5}=['13:45'];tmf{5}=['17:15'];

pathvl{6}='/oriole4/pu34/ac914/'
catchvl{6}='batch14.keep';tmn{6}=['17:15'];tmf{6}=['21'];

pathvl{7}='/oriole4/pu34/ac914/'
catchvl{7}='batch15.keep.rand';tmn{7}=['7'];tmf{7}=['8:57'];

pathvl{8}='/oriole4/pu34/lid915/'
catchvl{8}='batch.keep';tmn{8}=['8:57'];tmf{8}=['12:19'];

pathvl{9}='/oriole4/pu34/ac915/'
catchvl{9}='batch.keep';tmn{9}=['12:19'];tmf{9}=['15:32'];

pathvl{10}='/oriole4/pu34/500mu915/'
catchvl{10}='batch.keep';tmn{10}=['15:32'];tmf{10}=['18:19'];

pathvl{11}='/oriole4/pu34/ac915-2/'
catchvl{11}='batch15.keep';tmn{11}=['18:19'];tmf{11}=['21'];

pathvl{12}='/oriole4/pu34/ac915-2/'
catchvl{12}='batch16.keep.rand';tmn{12}=['7'];tmf{12}=['8:23'];

pathvl{13}='/oriole4/pu34/500mu916/'
catchvl{13}='batch.keep.rand';tmn{13}=['8:23'];tmf{13}=['11:50'];

pathvl{14}='/oriole4/pu34/ac916/'
catchvl{14}='batch.keep.rand';tmn{14}=['11:50'];tmf{14}=['21'];

pathvl{15}='/oriole4/pu34/wnon916/'
catchvl{15}='batch.catch.keep';tmn{15}=['7'];tmf{15}=['21'];

pathvl{16}='/oriole4/pu34/ac920/'
catchvl{16}='batch22.keep.rand';tmn{16}=['7'];tmf{16}=['11:29'];

pathvl{17}='/oriole4/pu34/mu922/'
catchvl{17}='batch22.keep.rand';tmn{17}=['11:29'];tmf{17}=['14:59'];

pathvl{18}='/oriole4/pu34/ac922/'
catchvl{18}='batch22.keep.rand';tmn{18}=['14:59'];tmf{18}=['21'];

pathvl{19}='/oriole4/pu34/wnon916/'
catchvl{19}='batch20.catch.keep';tmn{19}=['7'];tmf{19}=['11:51'];

pathvl{20}='/oriole4/pu34/500mu920/'
catchvl{20}='batch.keep.rand';tmn{20}=['11:51'];tmf{20}=['15:26'];

pathvl{21}='/oriole4/pu34/ac920/'
catchvl{21}='batch20.keep.rand';tmn{21}=['15:26'];tmf{21}=['21'];

pathvl{22}='/oriole4/pu34/ac922/'
catchvl{22}='batch24.keep.rand';tmn{22}=['7'];tmf{22}=['9:57'];

pathvl{23}='/oriole4/pu34/500mu924/'
catchvl{23}='batch.keep.rand';tmn{23}=['9:57'];tmf{23}=['13:51'];

pathvl{24}='/oriole4/pu34/ac924/'
catchvl{24}='batch.keep.rand';tmn{24}=['13:51'];tmf{24}=['21'];


pathvl{25}='/oriole4/pu34/ac924/'
catchvl{25}='batch25.keep.rand';tmn{25}=['7'];tmf{25}=['8:40'];

pathvl{26}='/oriole4/pu34/500mu925/'
catchvl{26}='batch.keep.rand';tmn{26}=['8:37'];tmf{26}=['12:42'];


pathvl{27}='/oriole4/pu34/iboacid925/'
catchvl{27}='batch';tmn{27}=['12:42'];tmf{27}=['21'];

pathvl{28}='/oriole4/pu34/ac925-2/'
catchvl{28}='batch26.keep';tmn{28}=['7'];tmf{28}=['10'];

pathvl{29}='/oriole4/pu34/ac926/'
catchvl{29}='batch05.keep.rand';tmn{29}=['7'];tmf{29}=['21'];

pathvl{30}='/oriole4/pu34/ac926/'
catchvl{30}='batch27.catch.keep.rand';tmn{30}=['7'];tmf{30}=['21'];

pathvl{31}='/oriole4/pu34/ac926/'
catchvl{31}='batch30.catch.keep';tmn{31}=['7'];tmf{31}=['21'];

pathvl{32}='/oriole4/pu34/ac926/'
catchvl{32}='batch02.catch.keep.rand';tmn{32}=['7'];tmf{32}=['21'];

pathvl{33}='/oriole4/pu34/screen1024/'
catchvl{33}='batchcomb';tmn{33}=['7'];tmf{33}=['21'];

pathvl{34}='/oriole4/pu34/ac926/'
catchvl{34}='batch29.catch.keep.rand';tmn{34}=['7'];tmf{34}=['21'];

pathvl{35}='/oriole4/pu34/ac926/'
catchvl{35}='batch28.catch.keep';tmn{35}=['7'];tmf{35}=['21'];


pathvl{36}='/oriole4/pu34/wnon1027/'
catchvl{36}='batchcomb27';tmn{36}=['7'];tmf{36}=['21'];

pathvl{37}='/oriole4/pu34/wnon1027/'
catchvl{37}='batchcomb31';tmn{37}=['7'];tmf{37}=['21'];

pathvl{38}='/oriole4/pu34/wnon1027/'
catchvl{38}='batchcomb02';tmn{38}=['7'];tmf{38}=['21'];


% pathvl{29}='/oriole4/pu34/dirsongs25/'
% catchvl{29}='batch';tmn{28}=['7'];tmf{28}=['21'];

ac(49)=[0];
avls.analind=[1:49]

%how to deal with dirfiles -- put all in one directory.
dirpathvl{1}='/oriole4/pu34/screen1024/'
dircatchvl{1}='batch.keepdir'
% 
% dirpathvl{2}='/oriole4/pk32bk28/wnon428/'
% dircatchvl{2}='batchcomb'
% 

numnotes=2
notes='ab'

fbins={};
tbinshft{1}=0.075;
tbinshft{2}=0.009;  
tbinshft{3}=0.009;
avls.pitchalign{1}=[0]
avls.pitchalign{2}=0;
avls.pitchalign{3}=0;



NFFT(1)=1024
NFFT(2)=256
NFFT(3)=256; 

fbins{1}=[6400,7400];
fbins{2}=[2000,2600];
fbins{3}=[2000,2600];

clear NT    
NT{1}='a'
NT{2}='b' 

 NT{3}='b'
avls.exclude=[0 1 1]
avls.pitchalign{1}=0;
avls.pitchalign{2}=1;
avls.pitchalighn{3}=1;
PRENT{1}='';PSTNT{1}='';
PRENT{2}='a';PSTNT{2}='b';
PRENT{3}='b';PSTNT{3}='b';
  %%%plotvals
  
acon=[ 1 3 4 6 7 9 11 12 14 15  16 18 19 21 22 24  25 27 28 29 30  31 32 33 34 35 36 37 38] ;  
muon=[2 5 8 10 13 17 20 23 26 ]
revanal{1}{1}=[];
muanal{1}{1}=[2 5 8 10 13 17 20 23] 




diron=[];
extraon=[]

clear colvals
colvals{1}=muon;
colvals{2}=acon
% colvals{3}=diron;

graphvals.numcol=2
graphvals.col='rk'
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
% avls.usex=usex
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
avls.mkfv=[38]
avls.PRENT=PRENT
avls.PSTNT=PSTNT
avls.repeatanal=[0 0 0]
avls.deadtm=.0833/24;
avls.muoffset=.9167/24;
avls.acoffset=.9167/24;
avls.pretm=4/24;


strcmd=['cd ' sumpath 'datasum']
eval(strcmd)

strcmd=['save ' matfilename '.mat avls graphvals'];
eval(strcmd);


