%diff with dat1 is that dat 2 uses wn files for mu trials.
%designed to call pitchsyntaxanal, and then inactivanal.m
clear date
clear avls
clear ac
clear graphvals
sumpath='/oriole5/pk16r47/'
strcmd=['cd ' sumpath]
eval(strcmd);
if (~exist('datasum','dir'))
    !mkdir datasum
end
matfilename='pathvals2'

extractvals=1


%white noise information
wnon{1}=  '2008-04-17  19:00:00'
minbnds{1}=5800
maxbnds{1}=6760
wnoff{1} ='2008-04-24 07:00:00'





clear pathvl vals




pathvl{1}='/oriole5/pk16r47/ac416/'
catchvl{1}='batch16.keep.rand'
timon=[]
timoff=[]
date{1}=''
% timon{1}='13:57:00'
% timoff{1}='18:51:00'
% date{1}='2008-01-26'

pathvl{2}='/oriole5/pk16r47/wnon417/'
catchvl{2}='batch20.catch.keep'


pathvl{3}='/oriole5/pk16r47/mu420/'
catchvl{3}='batchcomb'
timon{3}='14:47:00'
timoff{3}='18:56:00'
date{3}='2008-04-20'

pathvl{4}='/oriole5/pk16r47/ac420/'
catchvl{4}='batchcomb'

pathvl{5}='/oriole5/pk16r47/ac420/'
catchvl{5}='batch22.catch.keep'


pathvl{6}='/oriole5/pk16r47/300mu422/'
catchvl{6}='batchcomb'
timon{6}='14:47:00'
timoff{6}='18:56:00'
date{6}='2008-04-22'

pathvl{7}='/oriole5/pk16r47/ac422/'
catchvl{7}='batchcomb'

pathvl{8}='/oriole5/pk16r47/dirfiles/'
catchvl{8}='batch.catch'
timon{8}=''
timoff{8}=''
date{8}=''

pathvl{9}='/oriole5/pk16r47/dirfiles/'
catchvl{9}='batch'
timon{9}=''
timoff{9}=''
date{9}=''

pathvl{10}='/oriole5/pk16r47/wnon417/'
catchvl{10}='batch19.catch.keep'

pathvl{11}='/oriole5/pk16r47/250mu417/'
catchvl{11}='batch2'
timon{11}='11:15:00'
timoff{11}='15:15:00'
date{11}='2008-04-17'

pathvl{12}='/oriole5/pk16r47/ac417/'
catchvl{12}='batch.keep.rand'



pathvl{13}='/oriole5/pk16r47/lid416/'
catchvl{13}='batchcomb'
timon{13}='13:28:00'
timoff{13}='16:46:00'
date{13}='2008-04-16'

pathvl{14}='/oriole5/pk16r47/probein/'
catchvl{14}='batch16.keep'

pathvl{15}='/oriole5/pk16r47/ac416/'
catchvl{15}='batch17.keep'

timon{15}=''
timoff{15}=''
date{15}=''
usex(1:15)=0;
ac(15)=[0];
avls.analind=[1:15]

numnotes=2
notes='a'

fbins={};
tbinshft{1}=0.02;
  tbinshft{2}=0.02
%  tbinshft{3}=0.02
%  tbinshft{4}=.02
NFFT(1)=512
 NFFT(2)=512
% NFFT(3)=512
% NFFT(4)=512
fbins{1}=[6100,7000];
fbins{2}=[6100,7000];
% fbins{3}=[3000,4000];
% fbins{4}=[3000,4000];

clear NT    
NT{1}='e'
  NT{2}='z'
%  NT{3}='a'
%  NT{4}='a'
% NT{3}='g'

PRENT{1}='';PSTNT{1}='';
PRENT{2}='';PSTNT{2}='';
% PRENT{3}='b';PSTNT{3}='';
% PRENT{4}='';PSTNT{4}='';
  %%%plotvals
acon=[ 1 2 4 5 7 10 12 14 15];  
muon=[3 6 11 13 ]
diron=[8 9]
extraon=[]

clear colvals
colvals{1}=muon;
colvals{2}=acon
colvals{3}=diron;

graphvals.numcol=3
graphvals.col='rkc'
graphvals.plttext=1
graphvals.txtht=[3260 3200]
graphvals.tickht=8000
graphvals.edges{1}=[6000:60:7000]
graphvals.edges{2}=[6000:60:7000]

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
avls.mkfv=[13:15]
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


