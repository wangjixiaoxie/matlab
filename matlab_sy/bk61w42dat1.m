%designed to call pitchsyntaxanal, and then inactivanal.m
clear date
clear avls
clear ac tbinshft PRENT PSTNT muanal changenote
clear graphvals edges muanal revanal wn
clear tmn;clear tmf;
sumpath='/oriole6/bk61w42/'
strcmd=['cd ' sumpath]
eval(strcmd);
if (~exist('datasum','dir'))
    !mkdir datasum
end
matfilename='pathvals1'
clear pathvl vals wn;

%this is calculated with two days of canin data.
avls.initmean{1}=7022
avls.maxmeans=[6614 7139]
avls.initsd{1}=121
avls.exsongpath=1;
avls.songind=30;
avls.songtm=[2.05 3.05]

avls.initmean{2}=1703;
avls.initsd{2}=49.8;

avls.initmean{3}=1614;
avls.initsd{3}=40.6;

%%%TKTK
edges{1}=[6200:50:7900]
edges{2}=[1400:30:1900]
edges{3}=[1400:30:1900]

wn(1).tmon{1}={'2008-12-5 7'  }
wn(1).tmon{2}={'2008-12-20 7'}
wn(1).tmlist{1}={'2008-12-5 7' '2008-12-5 14:30' '2008-12-5 18:00' '2008-12-6 18:20' };
wn(1).ptlist{1}=[7115 7200 7350 7450]

wn(1).tmoff{1}={'2008-12-16 7'}
wn(1).tmoff{2}={'2008-12-25 7'}

wnrev(1).tmon={};
wnrev(1).tmoff={};
wnrev(2)=wnrev(1);

pathvl{1}='/oriole6/bk61w42/probin/'
catchvl{1}='batch30.keep';tmn{1}=['7'];tmf{1}=['09:22'];

pathvl{2}='/oriole6/bk61w42/lid1130/'
catchvl{2}='batch.keep';tmn{2}=['09:22'];tmf{2}=['13:22'];

pathvl{3}='/oriole6/bk61w42/ac1130/'
catchvl{3}='batch30.keep.rand';tmn{3}=['12:30'];tmf{3}=['21'];

pathvl{4}='/oriole6/bk61w42/ac1130/'
catchvl{4}='batchcomb';tmn{4}=['7'];tmf{4}=['10:24'];

pathvl{5}='/oriole6/bk61w42/lid1201/'
catchvl{5}='batchcomb';tmn{5}=['10:24'];tmf{5}=['14:22'];

pathvl{6}='/oriole6/bk61w42/ac1201/'
catchvl{6}='batchcomb';tmn{6}=['14:22'];tmf{6}=['21'];

pathvl{7}='/oriole6/bk61w42/ac1201/'
catchvl{7}='batch02.keep.rand';tmn{7}=['7'];tmf{7}=['9:13'];

pathvl{8}='/oriole6/bk61w42/500mu1202/'
catchvl{8}='batch';tmn{8}=['9:13'];tmf{8}=['12:56'];

pathvl{9}='/oriole6/bk61w42/ac1202/'
catchvl{9}='batch02.keep.rand';tmn{9}=['12:56'];tmf{9}=['21'];

pathvl{10}='/oriole6/bk61w42/ac1202/'
catchvl{10}='batch03.keep.rand';tmn{10}=['6'];tmf{10}=['8:10'];

pathvl{11}='/oriole6/bk61w42/500mu1203/'
catchvl{11}='batch.keep';tmn{11}=['8:13'];tmf{11}=['12:08'];

pathvl{12}='/oriole6/bk61w42/ac1203/'
catchvl{12}='batch03comb';tmn{12}=['12:08'];tmf{12}=['21'];

pathvl{13}='/oriole6/bk61w42/ac1203/'
catchvl{13}='batch04.keep';tmn{13}=['7'];tmf{13}=['8:27'];

pathvl{14}='/oriole6/bk61w42/400mu1204/'
catchvl{14}='batch.keep';tmn{14}=['8:27'];tmf{14}=['12:10'];

pathvl{15}='/oriole6/bk61w42/ac1204/'
catchvl{15}='batch04comb';tmn{15}=['12:10'];tmf{15}=['21'];

pathvl{16}='/oriole6/bk61w42/wnon1205newtemp/'
catchvl{16}='batch07comb';tmn{16}=['7'];tmf{16}=['9:45'];

pathvl{17}='/oriole6/bk61w42/400mu1207/'
catchvl{17}='batch.keep';tmn{17}=['9:45'];tmf{17}=['13:26'];

pathvl{18}='/oriole6/bk61w42/ac1207/'
catchvl{18}='batch07comb';tmn{18}=['13:26'];tmf{18}=['21'];


pathvl{19}='/oriole6/bk61w42/ac1207/'
catchvl{19}='batch09.keep.rand';tmn{19}=['7'];tmf{19}=['9:29'];

pathvl{20}='/oriole6/bk61w42/400mu1209/'
catchvl{20}='batch.keep';tmn{20}=['9:29'];tmf{20}=['13:29'];

pathvl{21}='/oriole6/bk61w42/ac1209/'
catchvl{21}='batch09';tmn{21}=['13:29'];tmf{21}=['21'];


pathvl{22}='/oriole6/bk61w42/ac1209/'
catchvl{22}='batch11.keep.rand';tmn{22}=['7'];tmf{22}=['8:50'];

pathvl{23}='/oriole6/bk61w42/400mu1211/'
catchvl{23}='batch.keep';tmn{23}='8:50';tmf{23}=['13:29'];

pathvl{24}='/oriole6/bk61w42/ac1211/'
catchvl{24}='batch11.keep';tmn{24}=['13:29'];tmf{24}=['21'];

pathvl{25}='/oriole6/bk61w42/ac1211/'
catchvl{25}='batch12.keep.rand';tmn{25}=['7'];tmf{25}=['9'];

pathvl{26}='/oriole6/bk61w42/400mu1212/'
catchvl{26}='batch';tmn{26}='9';tmf{26}=['12:29'];

pathvl{27}='/oriole6/bk61w42/ac1212/'
catchvl{27}='batch12.keep.rand';tmn{27}=['13:29'];tmf{27}=['21'];

pathvl{28}='/oriole6/bk61w42/ac1214/'
catchvl{28}='batch15.keep';tmn{28}=['7'];tmf{28}=['8:58'];

pathvl{29}='/oriole6/bk61w42/400mu1215/'
catchvl{29}='batch';tmn{29}=['8:58'];tmf{29}=['11:58'];

pathvl{30}='/oriole6/bk61w42/ac1215/'
catchvl{30}='batch';tmn{30}=['11:58'];tmf{30}=['21'];


pathvl{31}='/oriole6/bk61w42/ac1216/'
catchvl{31}='batch17.keep';tmn{31}=['7'];tmf{31}=['21'];

pathvl{32}='/oriole6/bk61w42/wnon1220/'
catchvl{32}='batch22.keep';tmn{32}=['6'];tmf{32}=['8:16'];

pathvl{33}='/oriole6/bk61w42/150mu1222/'
catchvl{33}='batch.keep';tmn{33}=['8:16'];tmf{33}=['11:28'];

pathvl{34}='/oriole6/bk61w42/ac1222/'
catchvl{34}='batch22.keep';tmn{34}=['11:29'];tmf{34}=['21'];

pathvl{35}='/oriole6/bk61w42/ac1222/'
catchvl{35}='batch23.keep';tmn{35}=['6'];tmf{35}=['7:28'];

pathvl{36}='/oriole6/bk61w42/150mu1223/'
catchvl{36}='batchcomb';tmn{36}=['7:28'];tmf{36}=['10:46'];

pathvl{37}='/oriole6/bk61w42/ac1223/'
catchvl{37}='batchcomb';tmn{37}=['10:46'];tmf{37}=['21'];

pathvl{38}='/oriole6/bk61w42/ac1218-2/'
catchvl{38}='batch19.keep';tmn{38}=['7'];tmf{38}=['7:38'];

pathvl{39}='/oriole6/bk61w42/150mu1219/'
catchvl{39}='batch';tmn{39}=['7:38'];tmf{39}=['10:26'];

pathvl{40}='/oriole6/bk61w42/ac1219/'
catchvl{40}='batch.catch';tmn{40}=['10:36'];tmf{40}=['21'];

pathvl{41}='/oriole6/bk61w42/wnon1205newtemp/'
catchvl{41}='batch05.catch.keep';tmn{41}=['7'];tmf{41}=['21'];

pathvl{42}='/oriole6/bk61w42/wnon1205newtemp/'
catchvl{42}='batch06comb.catch.keep';tmn{42}=['7'];tmf{42}=['21'];

pathvl{43}='/oriole6/bk61w42/ac1209/'
catchvl{43}='batch10';tmn{43}=['7'];tmf{43}=['21'];

pathvl{44}='/oriole6/bk61w42/ac1207/'
catchvl{44}='batch08';tmn{44}=['7'];tmf{44}=['21'];

pathvl{45}='/oriole6/bk61w42/ac1212/'
catchvl{45}='batch13';tmn{45}=['7'];tmf{45}=['21'];

pathvl{46}='/oriole6/bk61w42/wnon1220/'
catchvl{46}='batch21';tmn{46}=['7'];tmf{46}=['21'];

pathvl{47}='/oriole6/bk61w42/wnon1220/'
catchvl{47}='batch20.catch';tmn{47}=['7'];tmf{47}=['21'];


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
tbinshft{1}=0.04;
tbinshft{2}=0.024;  
tbinshft{3}=0.03;
avls.pitchalign{1}=[0]
avls.pitchalign{2}=0;
avls.pitchalign{3}=0;



NFFT(1)=512
NFFT(2)=512
NFFT(3)=512; 

fbins{1}=[5500 8000];
fbins{2}=[1200,1800];
fbins{3}=[1200,2000];


clear NT    
NT{1}='a'
NT{2}='b' 

 NT{3}='c'
avls.exclude=[0 1 1]
PRENT{1}='';PSTNT{1}='';
PRENT{2}='a';PSTNT{2}='';
PRENT{3}='';PSTNT{3}='';
  %%%plotvals
  
acon=[ 7 9 10 12 13 15 16 18 19 21 22 24 25 27 28 30 31 32  34 35 37 38 40:47] ;  
muon=[ 8 11 14 17 20 23 26 29 33 36 39]
revanal{1}{1}=[];
muanal{1}{1}=[ 8 11 14 17 20 23 26 29 33 36 39]

avls.changenote=1
avls.changeruns=[32:37]
avls.changebins=[5500 7200]



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
graphvals.wnrev=wnrev;
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
avls.mkfv=[40]
avls.PRENT=PRENT
avls.PSTNT=PSTNT
avls.repeatanal=[0 0 0]
avls.deadtm=.0833/24;
avls.muoffset=1/24;
avls.acoffset=1/24;
avls.pretm=4/24;


strcmd=['cd ' sumpath 'datasum']
eval(strcmd)

strcmd=['save ' matfilename '.mat avls graphvals'];
eval(strcmd);


