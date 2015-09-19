%designed to call pitchsyntaxanal, and then inactivanal.m
clear date
clear avls
clear ac tbinshft PRENT PSTNT muanal
clear graphvals edges muanal revanal wn
clear tmn;clear tmf;
sumpath='/oriole6/bk63w43/'
strcmd=['cd ' sumpath]
eval(strcmd);
if (~exist('datasum','dir'))
    !mkdir datasum
end
matfilename='pathvals1'
clear pathvl vals wn;

%this is calculated with two days of canin data.
avls.initmean{1}=6970
avls.maxmeans=[6614 7139]
avls.initsd{1}=128.1
avls.exsongpath=1;
avls.songind=30;
avls.songtm=[2.05 3.05]

avls.initmean{2}=1407;
avls.initsd{2}=60.8;

avls.initmean{3}=1614;
avls.initsd{3}=71.3;

%%%TKTK
edges{1}=[6300:55:7500]
edges{2}=[1800:60:2500]
edges{3}=[1800:60:2500]

wn(1).tmon{1}={'2008-12-5 7'}

wn(1).tmoff{1}={'2008-12-14 7'}

pathvl{1}='/oriole6/bk63w43/probin/'
catchvl{1}='batch01.keep.rand';tmn{1}=['7'];tmf{1}=['09:22'];

pathvl{2}='/oriole6/bk63w43/lid1201/'
catchvl{2}='batch.keep.rand';tmn{2}=['09:22'];tmf{2}=['13:22'];

pathvl{3}='/oriole6/bk63w43/ac1201/'
catchvl{3}='batch01.keep.rand';tmn{3}=['11:53'];tmf{3}=['15:53'];

pathvl{4}='/oriole6/bk63w43/ac1201/'
catchvl{4}='batch02.keep.rand';tmn{4}=['7'];tmf{4}=['10:24'];

pathvl{5}='/oriole6/bk63w43/500mu1202/'
catchvl{5}='batch02ac.keep';tmn{5}=['10:06'];tmf{5}=['14:00'];

pathvl{6}='/oriole6/bk63w43/500mu1202/'
catchvl{6}='batch02mu.keep';tmn{6}=['14:00'];tmf{6}=['21'];

pathvl{7}='/oriole6/bk63w43/500mu1202/'
catchvl{7}='batch03.keep';tmn{7}=['7'];tmf{7}=['9'];


pathvl{8}='/oriole6/bk63w43/400mu1203/'
catchvl{8}='batch.keep';tmn{8}=['9'];tmf{8}=['13'];
pathvl{9}='/oriole6/bk63w43/ac1203/'
catchvl{9}='batch.keep';tmn{9}=['13'];tmf{9}=['21'];

pathvl{10}='/oriole6/bk63w43/ac1203/'
catchvl{10}='batch04.keep';tmn{10}=['7'];tmf{10}=['9:38'];

pathvl{11}='/oriole6/bk63w43/400mu1204/'
catchvl{11}='batch04.keep';tmn{11}=['9:38'];tmf{11}=['13:27'];

pathvl{12}='/oriole6/bk63w43/ac1204/'
catchvl{12}='batch04.keep.rand';tmn{12}=['13:27'];tmf{12}=['21'];

pathvl{13}='/oriole6/bk63w43/wnon1205newtemp/'
catchvl{13}='batch08.catch.keep';tmn{13}=['7'];tmf{13}=['11:10'];

pathvl{14}='/oriole6/bk63w43/400mu1208/'
catchvl{14}='batch08.keep';tmn{14}=['11:10'];tmf{14}=['14:38'];

pathvl{15}='/oriole6/bk63w43/ac1208/'
catchvl{15}='batch08';tmn{15}=['14:38'];tmf{15}=['21'];

pathvl{16}='/oriole6/bk63w43/ac1208/'
catchvl{16}='batch09.keep.rand';tmn{16}=['7'];tmf{16}=['9:55'];

pathvl{17}='/oriole6/bk63w43/400mu1209/'
catchvl{17}='batch.keep';tmn{17}=['9:55'];tmf{17}=['13:37'];

pathvl{18}='/oriole6/bk63w43/ac1209/'
catchvl{18}='batch09.keep.rand';tmn{18}=['13:37'];tmf{18}=['21'];

pathvl{19}='/oriole6/bk63w43/ac1209/'
catchvl{19}='batch10.keep.rand';tmn{19}=['7'];tmf{19}=['9:10'];


pathvl{20}='/oriole6/bk63w43/250mu1210/'
catchvl{20}='batch10.keep';tmn{20}=['9:10'];tmf{20}=['11:37'];

pathvl{21}='/oriole6/bk63w43/ac1210/'
catchvl{21}='batch';tmn{21}=['13:37'];tmf{21}=['21'];

pathvl{22}='/oriole6/bk63w43/wnoffnewprobe/'
catchvl{22}='batch';tmn{22}=['7'];tmf{22}=['21'];



usex(1:49)=0;
ac(49)=[0];
avls.analind=[1:49]

%how to deal with dirfiles -- put all in one directory.
dirpathvl{1}='/oriole4/pk32bk28/dirfiles/'
dircatchvl{1}='batch'

dirpathvl{2}='/oriole4/pk32bk28/wnon428/'
dircatchvl{2}='batchcomb'


numnotes=2
notes='ab'

fbins={};
tbinshft{1}=0.04;
tbinshft{2}=.018
% tbinshft{2}=0.02;  
tbinshft{3}=0.03;
% avls.pitchalign{1}=[0]
% avls.pitchalign{2}=0;
% avls.pitchalign{3}=0;



NFFT(1)=512
NFFT(2)=512
NFFT(3)=256; 

fbins{1}=[6000 7200];
fbins{2}=[1000,1700];
fbins{3}=[1000,1900];

clear NT    
NT{1}='a'
NT{2}='b' 

 NT{3}='c'
avls.exclude=[0 1 1]
PRENT{1}='';PSTNT{1}='';
PRENT{2}='';PSTNT{2}='';
PRENT{3}='';PSTNT{3}='';
  %%%plotvals
  
acon=[ 1 3 4 6 7 9 10 12 13 15 16 18 19 21 22]  ;  
muon=[2 5 8 11 14 17 20]
revanal{1}{1}=[];
muanal{1}{1}=[2 5 8 11 14 17 20]
 


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
avls.mkfv=[22]
avls.PRENT=PRENT
avls.PSTNT=PSTNT
avls.repeatanal=[0 0 0]
avls.deadtm=.0833/24;
avls.muoffset=.75/24;
avls.acoffset=1/24;
avls.pretm=4/24;


strcmd=['cd ' sumpath 'datasum']
eval(strcmd)

strcmd=['save ' matfilename '.mat avls graphvals'];
eval(strcmd);


