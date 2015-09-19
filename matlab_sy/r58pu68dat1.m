%designed to call pitchsyntaxanal, and then inactivanal.m

sumpath='/doya2/r58pu68/'
strcmd=['cd ' sumpath]
eval(strcmd);
if (~exist('datasum','dir'))
    !mkdir datasum
end
matfilename='pathvals1'

extractvals=1

clear pathvl vals
clear avls
clear tbinshft

pathvl{1}='/doya2/r58pu68/postimplant/'
catchvl{1}='batchcomb3'
pathvl{2}='/doya2/r58pu68/templtest1/'
catchvl{2}='batchcomb2.keep.rand'
pathvl{3}='/doya2/r58pu68/tubesin/'
catchvl{3}='batch18comb'
pathvl{4}='/doya2/r58pu68/tubesin/'
catchvl{4}='batch1920comb.keep.rand'


usex(1:4)=0;
avls.repeatanal(1:4)=0
avls.analind=[1:4]
numnotes=1
notes='a'

fbins={};
tbinshft{1}=0.08;
NFFT(1)=1024
fbins{1}=[6000,8000];
clear NT    
NT{1}='a'


PRENT{1}='';PSTNT{1}='';

  %%%plotvals
avls.synind=[]
avls.translist{1}='e'
avls.translist{2}='c'
avls.translist{3}='z'
avls.translist{4}='i'
avls.translist{5}='k'
avls.pitchind=[1]
colvals='rmkbg'
graphvals.numcol=5
graphvals.col=colvals
graphvals.plttext=1
graphvals.txtht=[3260 3200]
graphvals.colvals=colvals
graphvals.pltpnts=1

graphvals.colvals=colvals

avls.pvls=pathvl
avls.datfile=matfilename;
avls.cvl=catchvl
avls.sumpath=sumpath
avls.mtflnm=matfilename  
avls.supanal=1
avls.NT=NT
avls.NFFT=NFFT
avls.fbins=fbins
avls.tshft=tbinshft
avls.usex=usex
avls.numnotes=numnotes
avls.mkfv=[1:4]
avls.PRENT=PRENT
avls.PSTNT=PSTNT
avls.bnds{1}='2007-09-20 07:00:00'
avvls.bnds{2}='2007-11-20 07:00:00'

strcmd=['cd ' sumpath 'datasum']
eval(strcmd)

strcmd=['save ' matfilename '.mat avls graphvals'];
eval(strcmd);


