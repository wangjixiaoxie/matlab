%designed to call pitchsyntaxanal, and then inactivanal.m
clear date
clear avls
clear ac
clear graphvals
sumpath='/oriole4/pk39bk5/'
strcmd=['cd ' sumpath]
eval(strcmd);
if (~exist('datasum','dir'))
    !mkdir datasum
end
matfilename='pathvals1'

extractvals=1


%white noise information
wnon{1}=  '2008-04-13  07:00:00'
minbnds{1}=6700
maxbnds{1}=7320
wnoff{1} ='2008-04-12 07:00:00'

wnon{2}='2008-04-14 07:00:00'
wnoff{2}='2008-03-14 07:00:00'
minbnds{2}=6700
maxbnds{2}=7450

wnon{3}='2008-03-14 07:00:00'
wnoff{3}='2008-03-17 07:00:00'
minbnds{3}=6700
maxbnds{3}=7630


wnon{4}='2008-03-19 14:00:00'
wnoff{4}='2008-03-22 07:00:00'
minbnds{4}=7280 
maxbnds{4}=8000

wnon{5}='2008-03-22 7:00:00'
wnoff{5}='2008-03-26 07:00:00'
minbnds{5}=6980
maxbnds{5}=8000




clear pathvl vals




pathvl{1}='/oriole4/pk39bk5/probein/'
catchvl{1}='batch27.keep'
timon=[]
timoff=[]
date{1}=''
% timon{1}='13:57:00'
% timoff{1}='18:51:00'
% date{1}='2008-01-26'

pathvl{2}='/oriole4/pk39bk5/lid427/'
catchvl{2}='batch.keep'
timon{2}='12:04:00'
timoff{2}='14:52:00'
date{2}='2008-04-27'

pathvl{3}='/oriole4/pk39bk5/ac427/'
catchvl{3}='batch.keep'


pathvl{4}='/oriole4/pk39bk5/temptest427/'
catchvl{4}='batch.keep.rand'

pathvl{5}='/oriole4/pk39bk5/200mu427/'
catchvl{5}='batchcomb2'
timon{5}='11:00:00'
timoff{5}='15:22:00'
date{5}='2008-04-28'

pathvl{6}='/oriole4/pk39bk5/ac428/'
catchvl{6}='batchcomb'

pathvl{7}='/oriole4/pk39bk5/ac428/'
catchvl{7}='batch29.keep.rand'

pathvl{8}='/oriole4/pk39bk5/300mu428newprobe/'
catchvl{8}='batchcomb'
timon{8}='13:00:00'
timoff{8}='15:15:00'
date{8}='2008-04-29'

pathvl{9}='/oriole4/pk39bk5/ac429/'
catchvl{9}='batch29.keep'
pathvl{10}='/oriole4/pk39bk5/wnon430/'
catchvl{10}='batch30.catch.keep'
pathvl{11}='/oriole4/pk39bk5/wnon430/'
catchvl{11}='batch01.catch.keep'

pathvl{12}='/oriole4/pk39bk5/wnon430/'
catchvl{12}='batch02.catch.keep'

pathvl{13}='/oriole4/pk39bk5/lid502/'
catchvl{13}='batch.catch.keep'

pathvl{14}='/oriole4/pk39bk5/300mu503/'
catchvl{14}='batch.catch.keep'


pathvl{15}='/oriole4/pk39bk5/ac502-2/'
catchvl{15}='batch03.catch.keep'

timon{15}=''
timoff{15}=''
date{15}=''
usex(1:15)=0;
ac(15)=[0];
avls.analind=[1:15]

numnotes=1
notes='a'

fbins={};
tbinshft{1}=0.09;
  tbinshft{2}=0.02
%  tbinshft{3}=0.02
%  tbinshft{4}=.02
NFFT(1)=1024
 NFFT(2)=512
% NFFT(3)=512
% NFFT(4)=512
fbins{1}=[6000,8000];
fbins{2}=[1000 2000]
% fbins{2}=[4700,5800];
% fbins{3}=[3000,4000];
% fbins{4}=[3000,4000];

clear NT    
NT{1}='a'
NT{2}='b'
%  NT{2}='b'
%  NT{3}='a'
%  NT{4}='a'
% NT{3}='g'

PRENT{1}='';PSTNT{1}='';
PRENT{2}='';PSTNT{2}='';
% PRENT{3}='b';PSTNT{3}='';
% PRENT{4}='';PSTNT{4}='';
  %%%plotvals
acon=[ 1 3 4 6 7 9 10 11 12 15];  
muon=[2 5 8 13 14]
diron=[]
extraon=[]

clear colvals
colvals{1}=muon;
colvals{2}=acon

graphvals.numcol=2
graphvals.col='rk'
graphvals.plttext=1
graphvals.txtht=[3260 3200]
graphvals.tickht=8000
graphvals.edges{1}=[6000:50:8000]
graphvals.edges{2}=[1000:50:2000]

graphvals.colvals=colvals
graphvals.pltpnts=1

graphvals.timon=timon;
graphvals.timoff=timoff;
graphvals.acon=ac;
graphvals.date=date;

graphvals.wnon=wnon;
graphvals.wnoff=wnoff;
graphvals.minbnds=minbnds;
graphvals.maxbnds=maxbnds;


graphvals.colvals=colvals
graphvals.chunkdata=1

avls.pvls=pathvl
avls.datfile=matfilename;
avls.cvl=catchvl
avls.sumpath=sumpath
avls.mtflnm=matfilename  
avls.supanal=0
avls.NT=NT
avls.NFFT=NFFT
avls.fbins=fbins
avls.tshft=tbinshft
avls.usex=usex
avls.numnotes=numnotes
avls.mkfv=[15]
avls.PRENT=PRENT
avls.PSTNT=PSTNT
avls.repeatanal=[0 0 0 0]
avls.bnds{1}='2007-09-20 07:00:00'
avvls.bnds{2}='2007-11-20 07:00:00'
avls.offset=2;
avls.acoffset=1;

strcmd=['cd ' sumpath 'datasum']
eval(strcmd)

strcmd=['save ' matfilename '.mat avls graphvals'];
eval(strcmd);


