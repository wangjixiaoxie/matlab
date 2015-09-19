%designed to call pitchsyntaxanal, and then inactivanal.m
clear date
clear avls
clear ac tbinshft PRENT PSTNT muanal changenote
clear graphvals edges muanal revanal wn
clear tmn;clear tmf;
sumpath='/oriole5/r32pu82/'
strcmd=['cd ' sumpath]
eval(strcmd);
if (~exist('datasum','dir'))
    !mkdir datasum
end
matfilename='pathvals1'
clear pathvl vals wn;

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
edges{1}=[6000:70:7700]
edges{2}=[1400:30:1900]
edges{3}=[1400:30:1900]

wn(1).tmon{1}={'2008-12-5 7'}

wn(1).tmoff{1}={'2008-12-14 7'}




pathvl{1}='/oriole5/r32pu82/screen/'
catchvl{1}='batch03.keep.rand';tmn{1}=['7'];tmf{1}=['21'];

pathvl{2}='/oriole5/r32pu82/probin/'
catchvl{2}='batch11.keep';tmn{2}=['7'];tmf{2}=['21'];

pathvl{3}='/oriole5/r32pu82/canin/'
catchvl{3}='batch07.keep';tmn{3}=['7'];tmf{3}=['21'];




usex(1:3)=0;
ac(3)=[0];
avls.analind=[1:3]

%how to deal with dirfiles -- put all in one directory.
dirpathvl{1}='/oriole4/pk32bk28/dirfiles/'
dircatchvl{1}='batch'

dirpathvl{2}='/oriole4/pk32bk28/wnon428/'
dircatchvl{2}='batchcomb'


numnotes=1
notes='ab'

fbins={};
tbinshft{1}=0.06;
% tbinshft{2}=0.024;  
% tbinshft{3}=0.03;
avls.pitchalign{1}=[0]
% avls.pitchalign{2}=0;
% avls.pitchalign{3}=0;



NFFT(1)=512
% NFFT(2)=512
% NFFT(3)=512; 

fbins{1}=[2000 3000];
% fbins{2}=[1200,1800];
% fbins{3}=[1200,2000];


clear NT    
NT{1}='a'
% NT{2}='b' 
% 
%  NT{3}='c'
avls.exclude=[0 1 1]
PRENT{1}='';PSTNT{1}='';
% PRENT{2}='';PSTNT{2}='';
% PRENT{3}='';PSTNT{3}='';
  %%%plotvals
  
acon=[ 1 2 3 ] ;  
muon=[]
revanal{1}{1}=[];
muanal{1}{1}=[]

avls.changenote=0
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
avls.mkfv=[1:3]
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


