%designed to call pitchsyntaxanal, and then inactivanal.m

sumpath='/doyale2/twarren/bk68b82/'
strcmd=['cd ' sumpath]
eval(strcmd);
if (~exist('datasum','dir'))
    !mkdir datasum
end
matfilename='pathvals1'

extractvals=1
makefv=1
clear pathvl vals
pathvl{1}='/doyale2/twarren/bk68b82/screen1/'
catchvl{1}='batch26dir'
pathvl{2}='/doyale2/twarren/bk68b82/screen1/'
catchvl{2}='batch29dir'
pathvl{3}='/doyale2/twarren/bk68b82/screen1/'
pathvl{4}='/doyale2/twarren/bk68b82/screen1/'

catchvl{3}='batch26.keep.rand'
catchvl{4}='batch29.keep.rand'

pathvl{5}='/doyale2/twarren/bk68b82/muon041007-1//'
catchvl{5}='batch.keep'
pathvl{6}='/doyale2/twarren/bk68b82/muon041007-2/'
catchvl{6}='batch.keep'
pathvl{7}='/doyale2/twarren/bk68b82/acon041007-1/'
catchvl{7}='batch.keep'
pathvl{8}='/doyale2/twarren/bk68b82/probeinacsf100907/'
catchvl{8}='batch.keep'
pathvl{9}='/doyale2/twarren/bk68b82/muon100907/'
catchvl{9}='batch.keep'
pathvl{10}='/doyale2/twarren/bk68b82/wna3/'
catchvl{10}='batch.catch.keep'


cd datasum
save pathvals.mat catchvl pathvl

usex(1:9)=0
numnotes=2
notes='ab'
clear fbins

tbinshft(1)=0.015;
    tbinshft(2)=0.015;
    NFFT(1)=1024
    NFFT(2)=1024;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
    fbins{1}=[2000,3000];
    fbins{2}=[3000 4000]
     NT{1}='a';PRENT{1}='';PSTNT{1}='';
    NT{2}='b';PRENT{2}='-';PSTNT{2}='b';
    
   
  %%%plotvals
  
muon=[5 6 9]
acon=[3 4 7 8 10]
diron=[1 2]
colvals{1}=diron;
colvals{2}=muon
colvals{3}=acon;

graphvals.numcol=3
graphvals.col='rmk'
graphvals.plttext=1
graphvals.txtht=[3260 3200]

supanal=0

 
    
strcmd=['save -append ' matfilename '.mat'];
    eval(strcmd);


