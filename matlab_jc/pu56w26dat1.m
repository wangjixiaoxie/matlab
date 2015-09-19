 %designed to call pitchsyntaxanal, and then inactivanal.m
clear date
clear avls
clear ac tbinshft PRENT PSTNT muanal
clear graphvals edges muanal revanal wn
clear tmn;clear tmf;
sumpath='/oriole5/pu56w26/'
strcmd=['cd ' sumpath]
eval(strcmd);
if (~exist('datasum','dir'))
    !mkdir datasum
end
matfilename='pathvals1'
clear pathvl vals wn;

%this is calculated with two days of canin data.
avls.initmean{1}=7056
avls.maxmeans=[6614 7139]
avls.initsd{1}=187
avls.exsongpath=1;
avls.songind=30;
avls.songtm=[2.05 3.05]

avls.initmean{2}=6657;
avls.initsd{2}=166;

% avls.initmean{3}=2220;
% avls.initsd{3}=168;

avls.changenote=0;
%%%TKTK
edges{1}=[5950:50:8500]
edges{2}=[1800:60:2500]
edges{3}=[1800:60:2500]

wn(1).tmon{1}={'2009-2-22 7' '2009-2-25 7'}
% wn(1).tmon{2}={'2008-7-27 7' '2008-7-30 7'}  
% wn(1).tmon{3}={'2008-9-19 7' }
wn(1).tmoff{1}={'2009-2-22 7' '2009-3-03 7'}
% wn(1).tmoff{2}={'2008-7-30 7' '2008-8-3 7'} 
% wn(1).tmoff{3}={'2008-9-26 7'}
wn(2).tmon=wn(1).tmon;
wn(2).tmoff=wn(1).tmoff;

% wn(3).tmon=wn(1).tmon;
% wn(3).tmoff=wn(1).tmoff;

wnrev(1).tmon{1}={'2009-3-03 7'}
wnrev(1).tmoff{1}={'2009-3-10 7'}
% wnrev(1).tmon{2}={'2009-3-3 7'}
% wnrev(1).tmoff{2}={'2009-3-10 7'}
wnrev(2).tmon=wnrev(1).tmon;
wnrev(2).tmoff=wnrev(1).tmoff;
% wnrev(3).tmon=wnrev(2).tmon;
% wnrev(3).tmoff=wnrev(2).tmoff;
graphvals.wnrev=wnrev;


pathvl{1}='/oriole5/pu56w26/21709-2/'
catchvl{1}='batch19.keep.rand';tmn{1}=['7'];tmf{1}=['9:07'];

pathvl{2}='/oriole5/pu56w26/21909_APV_2mM/'
catchvl{2}='batch.keep.rand';tmn{2}=['9:07'];tmf{2}=['12:42'];

pathvl{3}='/oriole5/pu56w26/219_ACSF/'
catchvl{3}='batch19';tmn{3}=['12:42'];tmf{3}=['21'];

pathvl{4}='/oriole5/pu56w26/224tmp3/'
catchvl{4}='batch25.keep.rand';tmn{4}=['7'];tmf{4}=['9:57'];

pathvl{5}='/oriole5/pu56w26/225_APV_2mM/'
catchvl{5}='batch25.keep.rand';tmn{5}=['9:57'];tmf{5}=['13:23'];

pathvl{6}='/oriole5/pu56w26/225ACSF/'
catchvl{6}='batch25.keep.rand';tmn{6}=['13:23'];tmf{6}=['21:00'];

pathvl{7}='/oriole5/pu56w26/225ACSF/'
catchvl{7}='batch27.keep.rand';tmn{7}=['7'];tmf{7}=['9:50'];

pathvl{8}='/oriole5/pu56w26/227_APV_2mM/'
catchvl{8}='batch27.keep';tmn{8}=['9:50'];tmf{8}=['13:18'];

pathvl{9}='/oriole5/pu56w26/227_ACSF/'
catchvl{9}='batch27.keep';tmn{9}=['13:18'];tmf{9}=['21'];


pathvl{10}='/oriole5/pu56w26/21709_ACSF/'
catchvl{10}='batch.keep';tmn{10}=['12'];tmf{10}=['15:55'];

pathvl{11}='/oriole5/pu56w26/21709_APV_2mM/'
catchvl{11}='batch.keep';tmn{11}=['15:55'];tmf{11}=['19:13'];

pathvl{12}='/oriole5/pu56w26/21709-2/'
catchvl{12}='batch17';tmn{12}=['19:13'];tmf{12}=['21'];


pathvl{13}='/oriole5/pu56w26/probin2_030109/'
catchvl{13}='batch02.keep';tmn{13}=['7'];tmf{13}=['10:10'];

pathvl{14}='/oriole5/pu56w26/0302_APV_2mM/'
catchvl{14}='batch.keep.rand';tmn{14}=['10:10'];tmf{14}=['13:15'];

pathvl{15}='/oriole5/pu56w26/030209_ACSF/'
catchvl{15}='batch.keep';tmn{15}=['13:15'];tmf{15}=['21'];

pathvl{16}='/oriole5/pu56w26/0303_wnrev/'
catchvl{16}='batch04.keep.rand';tmn{16}=['7'];tmf{16}=['10:24'];
pathvl{17}='/oriole5/pu56w26/0303_APV_2mM/'
catchvl{17}='batch.keep.rand';tmn{17}=['10:25'];tmf{17}=['13:50'];

pathvl{18}='/oriole5/pu56w26/0304_ACSF/'
catchvl{18}='batch04.keep.rand';tmn{18}=['13:50'];tmf{18}=['21'];


pathvl{19}='/oriole5/pu56w26/227_ACSF/'
catchvl{19}='batch01';tmn{19}=['7'];tmf{19}=['21'];

pathvl{20}='/oriole5/pu56w26/225ACSF/'
catchvl{20}='batch26';tmn{20}=['7'];tmf{20}=['21'];

pathvl{21}='/oriole5/pu56w26/224tmp3/'
catchvl{21}='batch24';tmn{21}=['7'];tmf{21}=['21'];

pathvl{22}='/oriole5/pu56w26/219_ACSF/'
catchvl{22}='batch20.catch.keep';tmn{22}=['7'];tmf{22}=['21'];

pathvl{23}='/oriole5/pu56w26/219_ACSF/'
catchvl{23}='batch21.catch.keep';tmn{23}=['7'];tmf{23}=['21'];

pathvl{24}='/oriole5/pu56w26/221ampon/'
catchvl{24}='batch22.catch.keep';tmn{24}=['7'];tmf{24}=['21'];

pathvl{25}='/oriole5/pu56w26/222tmp2/'
catchvl{25}='batch23.catch.keep';tmn{25}=['7'];tmf{25}=['21'];

usex(1:9)=0;
ac(49)=[0];
avls.analind=[1:9]

%how to deal with dirfiles -- put all in one directory.
dirpathvl{1}='/oriole4/pk32bk28/dirfiles/'
dircatchvl{1}='batch'

dirpathvl{2}='/oriole4/pk32bk28/wnon428/'
dircatchvl{2}='batchcomb'


numnotes=2
notes='ac'

fbins={};
tbinshft{1}=0.01;
tbinshft{2}=0.01;  
% tbinshft{3}=0.02;
avls.pitchalign{1}=[0]
avls.pitchalign{2}=0;
% avls.pitchalign{3}=0;



NFFT(1)=512
NFFT(2)=512
NFFT(3)=256; 

fbins{1}=[6000,7800];
fbins{2}=[6000,8400];
% fbins{3}=[6400,8500];

clear NT    
NT{1}='a'
NT{2}='c' 

%  NT{3}='c'
avls.exclude=[0 0 0]
PRENT{1}='';PSTNT{1}='';
PRENT{2}='';PSTNT{2}='';
% PRENT{3}='';PSTNT{3}='';
  %%%plotvals
  
acon=[ 1 3 4 6 7 9 10 12 13 15 16 18 19:25] ;  
muon=[2 5 8 11 14 17]
revanal{1}{1}=[];
revanal{1}{2}=[]
revanal{2}{1}=[]
revanal{2}{2}=[]
muanal{1}{1}=[2 5 8 11 14]


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
avls.mkfv=[3 22:25 ]
avls.PRENT=PRENT
avls.PSTNT=PSTNT
avls.repeatanal=[0 0 0]
avls.deadtm=.0833/24;
avls.muoffset=1.5/24;
avls.acoffset=1.5/24;
avls.pretm=4/24;


strcmd=['cd ' sumpath 'datasum']
eval(strcmd)

strcmd=['save ' matfilename '.mat avls graphvals'];
eval(strcmd);


