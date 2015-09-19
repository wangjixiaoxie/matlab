%designed to call pitchsyntaxanal, and then inactivanal.m
clear date
clear avls
clear ac
clear graphvals
clear timoff
clear timon
sumpath='/oriole2/pk2bk37/'
strcmd=['cd ' sumpath]
eval(strcmd);
if (~exist('datasum','dir'))
    !mkdir datasum
end
matfilename='pathvals1'

extractvals=1


%white noise information
wnon{1}=  '2008-04-29  07:00:00'
minbnds{1}=6900
maxbnds{1}=7430
wnoff{1} ='2008-5-14 07:00:00'

wnon{2}='2008-03-12 07:00:00'
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
pathvl{1}='/oriole2/pk2bk37/ac519/'
catchvl{1}='batch20.keep.rand'

pathvl{2}='/oriole2/pk2bk37/500mu520/'
catchvl{2}='batch'
timon{2}='10:42:00'
timoff{2}='13:45:00'
date{2}='2008-05-20'
pathvl{3}='/oriole2/pk2bk37/ac520/'
catchvl{3}='batch20'
timon{3}='13:45:00'
timoff{3}='19:45:00'
date{3}=date{2}

pathvl{4}='/oriole2/pk2bk37/ac520tmptest/'
catchvl{4}='batch21.keep.rand'

pathvl{5}='/oriole4/pk2bk37/750mu521/'
catchvl{5}='batch.keep'
timon{5}='09:45:00'
timoff{5}='14:00:00'
date{5}='2008-05-21'

pathvl{6}='/oriole2/pk2bk37/ac521/'
catchvl{6}='batch.keep.rand'
timon{6}='14:15:00'
timoff{6}='18:30:00'
date{6}='2008-05-21'

pathvl{7}='/oriole4/pk2bk37/750mu522/'
catchvl{7}='batch'

pathvl{8}='/oriole2/pk2bk37/ac522newprobe/'
catchvl{8}='batch.keep'

pathvl{9}='/oriole4/pk2bk37/1500mu523/'
catchvl{9}='batch'




pathvl{10}='/oriole2/pk2bk37/524newprobe/'
catchvl{10}='batch'

pathvl{11}='/oriole2/pk2bk37/500mu524/'
catchvl{11}='batch'

pathvl{12}='/oriole4/pk2bk37/ac524/'
catchvl{12}='batch26.keep.rand'

pathvl{13}='/oriole2/pk2bk37/lid525/'
catchvl{13}='batch'

pathvl{14}='/oriole4/pk2bk37/ac524/'
catchvl{14}='batch25comb'

pathvl{15}='/oriole2/pk2bk37/1000mu526/'
catchvl{15}='batch'

pathvl{16}='/oriole2/pk2bk37/ac526-2/'
catchvl{16}='batch.keep.rand'

pathvl{17}='/oriole2/pk2bk37/ac526/'
catchvl{17}='batch.keep'

pathvl{18}='/oriole2/pk2bk37/ampon527/'
catchvl{18}='batch29.keep'

pathvl{19}='/oriole4/pk2bk37/1000mu529/'
catchvl{19}='batch.catch'

pathvl{20}='/oriole2/pk2bk37/ac529/'
catchvl{20}='batch.catch'

timon{20}=''
timoff{20}=''
date{20}=''
usex(1:20)=0;
ac(20)=[0];
avls.analind=[1:20]

numnotes=1
notes='a'

fbins={};
tbinshft{1}=0.08;
  

NFFT(1)=1024
 

fbins{1}=[6000,8000];


clear NT    
NT{1}='a'
 
% NT{3}='g'

PRENT{1}='';PSTNT{1}='';

  %%%plotvals
acon=[ 1 3 4 6 8 10 12 14 16 17 18 20] ;  
muon=[2 5 7 9 11 13 15 19]
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
graphvals.edges{1}=[6700:80:7900]
graphvals.edges{2}=[2000:40:2500]
graphvals.edges{3}=[1300:40:1800]
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
avls.mkfv=[17:20]
avls.PRENT=PRENT
avls.PSTNT=PSTNT
avls.repeatanal=[0 0 0]
avls.bnds{1}='2007-09-20 07:00:00'
avvls.bnds{2}='2007-11-20 07:00:00'
avls.offset=1.5;
avls.acoffset=1;

strcmd=['cd ' sumpath 'datasum']
eval(strcmd)

strcmd=['save ' matfilename '.mat avls graphvals'];
eval(strcmd);


