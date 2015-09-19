%designed to call pitchsyntaxanal, and then inactivanal.m
clear date
clear avls
clear ac
clear graphvals
sumpath='/oriole5/bk15bk14/'
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




pathvl{1}='/oriole5/bk15bk14/probein/'
catchvl{1}='batch13.keep'
timon=[]
timoff=[]
date{1}=''
% timon{1}='13:57:00'
% timoff{1}='18:51:00'
% date{1}='2008-01-26'

pathvl{2}='/oriole5/bk15bk14/lid41308/'
catchvl{2}='batch.keep'
timon{2}='09:55:00'
timoff{2}='12:47:00'
date{2}='2008-04-13'

pathvl{3}='/oriole5/bk15bk14/acsf41308/'
catchvl{3}='batch13.keep.rand'


pathvl{4}='/oriole5/bk15bk14/acsf41308/'
catchvl{4}='batch14.keep.rand'

pathvl{5}='/oriole5/bk15bk14/250mu41408/'
catchvl{5}='batchcomb2'
timon{5}='13:53:00'
timoff{5}='17:30:00'
date{5}='2008-03-08'

pathvl{6}='/oriole5/bk15bk14/acsf41408/'
catchvl{6}='batch14.keep'

pathvl{7}='/oriole5/bk15bk14/wnon417newtemp/'
catchvl{7}='batch18comb'

pathvl{8}='/oriole5/bk15bk14/lid417/'
catchvl{8}='batch.catch.keep'
timon{8}='14:35:00'
timoff{8}='17:38:00'
date{8}='2008-04-17'

pathvl{9}='/oriole5/bk15bk14/ac418/'
catchvl{9}='batch18comb'



pathvl{10}='/oriole5/bk15bk14/ac418/'
catchvl{10}='batch19.catch.keep'

pathvl{11}='/oriole5/bk15bk14/250mu41908/'
catchvl{11}='batch.catch.keep'
timon{11}='11:45:00'
timoff{11}='16:00:00'
date{11}='2008-04-19'

pathvl{12}='/oriole5/bk15bk14/ac41908/'
catchvl{12}='batch19.catch.keep'



timon{12}=''
timoff{12}=''
date{12}=''
usex(1:12)=0;
ac(12)=[0];
avls.analind=[1:12]

numnotes=1
notes='a'

fbins={};
tbinshft{1}=0.02;
  tbinshft{2}=0.02
%  tbinshft{3}=0.02
%  tbinshft{4}=.02
NFFT(1)=1024
 NFFT(2)=1024
% NFFT(3)=512
% NFFT(4)=512
fbins{1}=[4700,5800];
% fbins{2}=[4700,5800];
% fbins{3}=[3000,4000];
% fbins{4}=[3000,4000];

clear NT    
NT{1}='b'
%  NT{2}='b'
%  NT{3}='a'
%  NT{4}='a'
% NT{3}='g'

PRENT{1}='';PSTNT{1}='';
% PRENT{2}='b';PSTNT{2}='';
% PRENT{3}='b';PSTNT{3}='';
% PRENT{4}='';PSTNT{4}='';
  %%%plotvals
acon=[ 1 3 4 6 7 9 10 12];  
muon=[2 5 8 11 ]
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
graphvals.edges{1}=[4000:80:6000]
graphvals.edges{2}=[4000:80:6000]
graphvals.edges{3}=[3000:50:4000]
graphvals.edges{4}=[3000:50:4000]
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
avls.mkfv=[1:12]
avls.PRENT=PRENT
avls.PSTNT=PSTNT
avls.repeatanal=[0 0 0 0]
avls.bnds{1}='2007-09-20 07:00:00'
avvls.bnds{2}='2007-11-20 07:00:00'
avls.offset=1.75/24;
avls.acoffset=1/24;
avls.deadtm=.2/24

strcmd=['cd ' sumpath 'datasum']
eval(strcmd)

strcmd=['save ' matfilename '.mat avls graphvals'];
eval(strcmd);


