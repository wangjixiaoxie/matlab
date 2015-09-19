%designed to call pitchsyntaxanal, and then inactivanal.m

sumpath='/doya/bk100bk99/'
strcmd=['cd ' sumpath]
eval(strcmd);
if (~exist('datasum','dir'))
    !mkdir datasum
end
matfilename='pathvals1'

extractvals=1

clear pathvl vals
pathvl{1}='/doya/bk100bk99/probein/'
catchvl{1}='batch2'
pathvl{2}='/doya/bk100bk99/probein/'
catchvl{2}='batch.keep'

pathvl{3}='/doya/bk100bk99/mu200_010608/'
catchvl{3}='batch'
timon{3}='10:40:00'
timoff{3}='15:00:00'
date{3}='2008-01-06'

pathvl{4}='/doya/bk100bk99/ac010608/'
catchvl{4}='batch.keep.rand'
timon{4}='12:44:00'
timoff{4}='16:47:00'
date{4}='2008-01-07'
pathvl{5}='/doya/bk100bk99/mu300_010708/'
catchvl{5}='batch.keep'
pathvl{6}='/doya/bk100bk99/ac010708/'
catchvl{6}='batch.keep'
pathvl{7}='/doya/bk100bk99/ac010808-temptest/'
catchvl{7}='batch'
pathvl{8}='/doya/bk100bk99/mu400_010808/'
catchvl{8}='batch'
timon{8}='12:00:00'
timoff{8}='16:00:00'
date{8}='2008-01-08'
pathvl{9}='/doya/bk100bk99/ac010808-2/'
catchvl{9}='batch'




usex(1:9)=0;
avls.analind=[1:9]

numnotes=3
notes='ab'

fbins={};
tbinshft{1}=0.015;
tbinshft{2}=.01
tbinshft{3}=.01
NFFT(1)=512
NFFT(2)=1024
NFFT(3)=1024
fbins{1}=[3000,4000];
    
fbins{2}=[4000 6000]
fbins{3}=[4000 6000]

clear NT    
NT{1}='b'
NT{2}='a'
NT{3}='a'



PRENT{1}='';PSTNT{1}='b';
PRENT{2}='';PSTNT{2}='a'
PRENT{3}='';PSTNT{3}='b'

   
  %%%plotvals
muon=[3 5 8]
acon=[1 2 4 6 7 9]
diron=[]
extraon=[]


colvals{1}=muon
colvals{2}=acon;


graphvals.numcol=2
graphvals.col='rk'
graphvals.plttext=1
graphvals.txtht=[3260 3200]
graphvals.colvals=colvals
graphvals.pltpnts=1

graphvals.timon=timon;
graphvals.timoff=timoff;
graphvals.date=date;
graphvals.colvals=colvals
graphvals.chunkdata=0

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
avls.mkfv=[9]
avls.PRENT=PRENT
avls.PSTNT=PSTNT
avls.repeatanal=[1 0 0]
avls.bnds{1}='2007-09-20 07:00:00'
avvls.bnds{2}='2007-11-20 07:00:00'

strcmd=['cd ' sumpath 'datasum']
eval(strcmd)

strcmd=['save ' matfilename '.mat avls graphvals'];
eval(strcmd);


