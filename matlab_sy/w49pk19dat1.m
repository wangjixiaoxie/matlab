%designed to call pitchsyntaxanal, and then inactivanal.m

sumpath='/oriole2/w49pk19/'
strcmd=['cd ' sumpath]
eval(strcmd);
if (~exist('datasum','dir'))
    !mkdir datasum
end
matfilename='pathvals1'

extractvals=1


clear pathvl vals
pathvl{1}='/oriole2/w49pk19/ac3408temptest/'
catchvl{1}='batch05.keep.rand'
timon=[]
timoff=[]
date{1}=''
% timon{1}='13:57:00'
% timoff{1}='18:51:00'
% date{1}='2008-01-26'

pathvl{2}='/oriole2/w49pk19/200mu3.5.08/'
catchvl{2}='batch.keep.rand'

pathvl{3}='/oriole2/w49pk19/ac3.5.08/'
catchvl{3}='batch.keep.rand'

pathvl{4}='/oriole2/w49pk19/280mu3.6.08/'
catchvl{4}='batch.keep.rand'


pathvl{5}='/oriole2/w49pk19/280mu3.6.08/'
catchvl{5}='batch2.keep'

pathvl{6}='/oriole2/w49pk19/280mu3.6.08/'
catchvl{6}='batch3.keep.rand'


pathvl{7}='/oriole2/w49pk19/ac3.5.08/'
catchvl{7}='batch06.keep.rand'



pathvl{8}='/oriole2/w49pk19/ac3608/'
catchvl{8}='batch.keep.rand'

pathvl{9}='/oriole2/w49pk19/ac3608/'
catchvl{9}='batch2.keep.rand'

pathvl{10}='/oriole2/w49pk19/temptest/'
catchvl{10}='batch01.rand.keep'

pathvl{11}='/oriole2/w49pk19/ac3608/'
catchvl{11}='batch2.keep.rand'

pathvl{12}='/oriole2/w49pk19/ampon306/'
catchvl{12}='batch.catch'

pathvl{13}='/oriole2/w49pk19/300mu3808/'
catchvl{13}='batch.catch'

pathvl{14}='/oriole2/w49pk19/300mu3908/'
catchvl{14}='batch.catch'

pathvl{15}='/oriole2/w49pk19/ac3908/'
catchvl{15}='batch.catch'

pathvl{16}='/oriole2/w49pk19/ac3808/'
catchvl{16}='batch.catch.keep'


pathvl{17}='/oriole2/w49pk19/200mu31108/'
catchvl{17}='batch.catch'

pathvl{18}='/oriole2/w49pk19/ac3908/'
catchvl{18}='batch.catch.keep.rand'

pathvl{19}='/oriole2/w49pk19/ac31108newtemp/'
catchvl{19}='batch12.keep.rand'

pathvl{20}='/oriole2/w49pk19/200mu31208ampoff/'
catchvl{20}='batch.keep.rand'

pathvl{21}='/oriole2/w49pk19/lid31308/'
catchvl{21}='batch.keep.rand'

pathvl{22}='/oriole2/w49pk19/ac31308_2/'
catchvl{22}='batch.keep.rand'


pathvl{23}='/oriole2/w49pk19/ac31308_2/'
catchvl{23}='batch14.keep.rand'

pathvl{24}='/oriole2/w49pk19/lid31408/'
catchvl{24}='batch3'

pathvl{25}='/oriole2/w49pk19/ac31408wnon/'
catchvl{25}='batch.catch'

pathvl{26}='/oriole2/w49pk19/ac31408wnon/'
catchvl{26}='batch16.catch.keep'

pathvl{27}='/oriole2/w49pk19/lido31608/'
catchvl{27}='batch.catch'

pathvl{28}='/oriole2/w49pk19/ac31608/'
catchvl{28}='batch.catch'

pathvl{29}='/oriole2/w49pk19/lid31708/'
catchvl{29}='batch.catch'



date{29}=''
usex(1:29)=0;
avls.analind=[1:29]

numnotes=2
notes='a'

fbins={};
tbinshft{1}=0.025;
 tbinshft{2}=0.025
% tbinshft{3}=0.015
NFFT(1)=512
 NFFT(2)=512
% NFFT(3)=512

fbins{1}=[6000,8000];
fbins{2}=[6000,8000];
% fbins{3}=[3000,4000];

clear NT    
NT{1}='a'
NT{2}='a'
% NT{2}='f'
% NT{3}='g'

PRENT{1}='-';PSTNT{1}='a';
PRENT{2}='-a';PSTNT{2}='';
PRENT{3}='';PSTNT{3}='';
   
  %%%plotvals
muon=[2 4 5 6 13 14 17 20 21 24 27 29];  
acon=[1 3 7 8 9 10 11 12 15 16 18 19 22 23 25 26 28]
diron=[]
extraon=[]

clear colvals
colvals{1}=muon;
colvals{2}=acon

graphvals.numcol=2
graphvals.col='rk'
graphvals.plttext=1
graphvals.txtht=[3260 3200]
graphvals.edges=[2500 4000]
graphvals.colvals=colvals
graphvals.pltpnts=1

graphvals.timon=timon;
graphvals.timoff=timoff;
graphvals.date=date;
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
avls.mkfv=[28:29]
avls.PRENT=PRENT
avls.PSTNT=PSTNT
avls.repeatanal=[0 0 0]
avls.bnds{1}='2007-09-20 07:00:00'
avvls.bnds{2}='2007-11-20 07:00:00'
avls.offset=2;

strcmd=['cd ' sumpath 'datasum']
eval(strcmd)

strcmd=['save ' matfilename '.mat avls graphvals'];
eval(strcmd);


