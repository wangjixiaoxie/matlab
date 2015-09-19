%designed to call pitchsyntaxanal, and then inactivanal.m

sumpath='/doya/pk11r11/'
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

%11-8%
pathvl{1}='/doya/pk11r11/templtest2/'
catchvl{1}='batch.keep.rand'
%11-10%
pathvl{2}='/doyale4/twarren/pk11r11/wnon/'
catchvl{2}='batch16.keep.catch.rand'

pathvl{3}='/doyale4/twarren/pk11r11/wnon/'
catchvl{3}='batch18.keep.catch'
usex(1:3)=0;
avls.repeatanal(1:13)=0
avls.analind=[1:3]
numnotes=2
notes='de'

fbins={};
tbinshft{1}=0.1;
tbinshft{2}=.015
NFFT(1)=1024
NFFT(2)=1024
fbins{1}=[2500,3600];
fbins{2}=fbins{1}  
clear NT    
NT{1}='d'
NT{2}='e'
% NT{3}='e'

PRENT{1}='';PSTNT{1}='';
PRENT{2}='';PSTNT{2}='';
PRENT{3}='';PSTNT{3}='';
%%%plotvals
avls.synind=[1]
avls.translist{1}='e'
avls.translist{2}='d'
avls.translist{3}='i'
%this is for the number of notes
avls.pitchind=[2]
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
avls.mkfv=[3]
avls.PRENT=PRENT
avls.PSTNT=PSTNT
avls.bnds{1}='2007-09-20 07:00:00'
avvls.bnds{2}='2007-11-20 07:00:00'

strcmd=['cd ' sumpath 'datasum']
eval(strcmd)

strcmd=['save ' matfilename '.mat avls graphvals'];
eval(strcmd);


