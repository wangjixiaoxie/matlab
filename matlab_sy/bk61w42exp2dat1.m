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
matfilename='pathvals3'
clear pathvl vals wn colvals;

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
edges{1}=[6000:45:8000]
edges{2}=[1400:14:1900]
edges{3}=[1400:30:1900]

wn(1).tmon{1}={'2009-03-22 7' '2009-03-27 7' }
% wn(1).tmon{2}={'2008-12-20 7'}

wn(1).tmoff{1}={'2009-03-26 7' '2009-04-05 7'}
% wn(1).tmoff{2}={'2008-12-25 7'}

wnrev(1).tmon={};
wnrev(1).tmoff={};
wnrev(2)=wnrev(1);

pathvl{1}='/oriole6/bk61w42/wnoff320pm/'
catchvl{1}='batch21.keep';tmn{1}=['7'];tmf{1}=['10:17'];

pathvl{2}='/oriole6/bk61w42/125mu321/'
catchvl{2}='batch21.keep';tmn{2}=['10:17'];tmf{2}=['13:23'];

pathvl{3}='/oriole6/bk61w42/ac321/'
catchvl{3}='batch.keep';tmn{3}=['13:23'];tmf{3}=['21'];

pathvl{4}='/oriole6/bk61w42/wnon322/'
catchvl{4}='batch25.keep';tmn{4}=['7'];tmf{4}=['9:32'];

pathvl{5}='/oriole6/bk61w42/125mu325/'
catchvl{5}='batch.keep';tmn{5}=['9:32'];tmf{5}=['12:39'];

pathvl{6}='/oriole6/bk61w42/ac325/'
catchvl{6}='batch25.keep.rand';tmn{6}=['12:39'];tmf{6}=['21'];

pathvl{7}='/oriole6/bk61w42/ac325/'
catchvl{7}='batch27.keep';tmn{7}=['7'];tmf{7}=['9:57'];

pathvl{8}='/oriole6/bk61w42/125mu327/'
catchvl{8}='batch27.keep';tmn{8}=['9:57'];tmf{8}=['13:10'];

pathvl{9}='/oriole6/bk61w42/ac327/'
catchvl{9}='batch27.keep.rand';tmn{9}=['13:10'];tmf{9}=['21'];

pathvl{10}='/oriole6/bk61w42/ac327/'
catchvl{10}='batch30.keep.rand';tmn{10}=['7'];tmf{10}=['9:41'];

pathvl{11}='/oriole6/bk61w42/125mu330/'
catchvl{11}='batch30.keep';tmn{11}=['9:41'];tmf{11}=['12:57'];

pathvl{12}='/oriole6/bk61w42/ac330/'
catchvl{12}='batch';tmn{12}=['12:57'];tmf{12}=['21'];

pathvl{13}='/oriole6/bk61w42/ac330/'
catchvl{13}='batch01.keep';tmn{13}=['7'];tmf{13}=['9:44'];

pathvl{14}='/oriole6/bk61w42/125mu401/'
catchvl{14}='batch01.keep';tmn{14}=['9:44'];tmf{14}=['13:12'];

pathvl{15}='/oriole6/bk61w42/ac401/'
catchvl{15}='batch01.keep';tmn{15}=['12:57'];tmf{15}=['21'];

pathvl{16}='/oriole6/bk61w42/ac401/'
catchvl{16}='batch03.keep';tmn{16}=['7'];tmf{16}=['8:57'];

pathvl{17}='/oriole6/bk61w42/125mu403/'
catchvl{17}='batch.keep';tmn{17}=['8:57'];tmf{17}=['11:47'];

pathvl{18}='/oriole6/bk61w42/ac403-2/'
catchvl{18}='batch.keep.rand';tmn{18}=['11:47'];tmf{18}=['21'];


pathvl{19}='/oriole6/bk61w42/wnon322/'
catchvl{19}='batch23.catch.keep';tmn{19}=['7'];tmf{19}=['21'];


pathvl{20}='/oriole6/bk61w42/wnon322/'
catchvl{20}='batch24.catch.keep';tmn{20}=['7'];tmf{20}=['21'];


pathvl{21}='/oriole6/bk61w42/ac325/'
catchvl{21}='batch26.catch.keep';tmn{21}=['7'];tmf{21}=['21'];

pathvl{22}='/oriole6/bk61w42/ac327/'
catchvl{22}='batch29';tmn{22}=['7'];tmf{22}=['21'];

pathvl{23}='/oriole6/bk61w42/ac330/'
catchvl{23}='batch31';tmn{23}=['7'];tmf{23}=['21'];

pathvl{24}='/oriole6/bk61w42/ac401/'
catchvl{24}='batch02';tmn{24}=['7'];tmf{24}=['21'];


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

fbins{1}=[6500 8000];
fbins{2}=[1200,2200];
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
  
acon=[1 3 4 6 7 9 10 12 13 15 16 18:24] ;  
muon=[2 5 8 11 14 17]
revanal{1}{1}=[];
muanal{1}{1}=[ 2 5 8 11 14 17 ]

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
avls.mkfv=[22:24]
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


