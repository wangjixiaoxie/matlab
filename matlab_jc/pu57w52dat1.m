%designed to call pitchsyntaxanal, and then inactivanal.m
clear date
clear avls
clear ac tbinshft PRENT PSTNT muanal
clear graphvals edges muanal revanal wn
clear tmn;clear tmf;
sumpath='/oriole5/pu57w52/'
strcmd=['cd ' sumpath]
eval(strcmd);
if (~exist('datasum','dir'))
    !mkdir datasum
end
matfilename='pathvals1'
clear pathvl vals wn wnrev;

%this is calculated with two days of canin data.
avls.initmean{1}=7185
avls.maxmeans=[6614 7139]
avls.initsd{1}=133
avls.exsongpath=1;
avls.songind=30;
avls.songtm=[2.05 3.05]

avls.initmean{2}=6746;
avls.initsd{2}=113;

% avls.initmean{3}=2220;
% avls.initsd{3}=168;

avls.changenote=0;
%%%TKTK
edges{1}=[6000:50:8500]
edges{2}=[1800:60:2500]
edges{3}=[1800:60:2500]

wn(1).tmon{1}={'2009-2-10 7' '2009-2-12 7'}
% wn(1).tmon{2}={'2009-2-10 7' '2009-2-12 7'}  
% wn(1).tmon{3}={'2008-9-19 7' }
wn(1).tmoff{1}={'2009-2-12 7' '2009-2-18 7'}
% wn(1).tmoff{2}={'2009-2-12 7' '2009-2-18 7'} 
% wn(1).tmoff{3}={'2008-9-26 7'}
wn(2).tmon=wn(1).tmon;
wn(2).tmoff=wn(1).tmoff;
% 
% wn(3).tmon=wn(1).tmon;
% wn(3).tmoff=wn(1).tmoff;

wnrev(1).tmon{1}={'2009-2-18 7'}
wnrev(1).tmoff{1}={'2009-2-22 7'}
% wnrev(1).tmon{2}={'2009-2-18 7'}
% wnrev(1).tmoff{2}={'2009-2-22 7'}
wnrev(2).tmon=wnrev(1).tmon;
wnrev(2).tmoff=wnrev(1).tmoff;
% wnrev(3).tmon=wnrev(2).tmon;
% wnrev(3).tmoff=wnrev(2).tmoff;
graphvals.wnrev=wnrev;


pathvl{1}='/oriole5/pu57w52/20408_AC/'
catchvl{1}='batch06.keep.rand';tmn{1}=['7'];tmf{1}=['10'];

pathvl{2}='/oriole5/pu57w52/20609_APV_2mM/'
catchvl{2}='batch.keep.rand';tmn{2}=['10'];tmf{2}=['14:00'];

pathvl{3}='/oriole5/pu57w52/20709_AC/'
catchvl{3}='batch08.keep.rand';tmn{3}=['7'];tmf{3}=['9:41'];

pathvl{4}='/oriole5/pu57w52/20809_APV10mM/'
catchvl{4}='batch1.keep.rand';tmn{4}=['9:41'];tmf{4}=['13:00'];

pathvl{5}='/oriole5/pu57w52/20809_ACSF/'
catchvl{5}='batch08.keep.rand';tmn{5}=['13'];tmf{5}=['21:00'];

pathvl{6}='/oriole5/pu57w52/20609_AC/'
catchvl{6}='batch';tmn{6}=['14'];tmf{6}=['21:00'];

pathvl{7}='/oriole5/pu57w52/20709_AC/'
catchvl{7}='batch07.keep.rand';tmn{7}=['16'];tmf{7}=['21'];

pathvl{8}='/oriole5/pu57w52/20609_AC/'
catchvl{8}='batch07';tmn{8}=['7'];tmf{8}=['12:00'];

pathvl{9}='/oriole5/pu57w52/20709_APV_5mM/'
catchvl{9}='batch';tmn{9}=['12'];tmf{9}=['16'];

pathvl{10}='/oriole5/pu57w52/20809_ACSF/'
catchvl{10}='batch09.keep';tmn{10}=['7'];tmf{10}=['9:25'];

pathvl{11}='/oriole5/pu57w52/20909_APV_5mM/'
catchvl{11}='batch.keep.rand';tmn{11}=['9:25'];tmf{11}=['13:00'];

pathvl{12}='/oriole5/pu57w52/20909_AC/'
catchvl{12}='batch.keep.rand';tmn{12}=['13'];tmf{12}=['21'];

pathvl{13}='/oriole5/pu57w52/wnon2909eve/'
catchvl{13}='batch12.keep';tmn{13}=['7'];tmf{13}=['11:44'];

pathvl{14}='/oriole5/pu57w52/21209_APV_5mM/'
catchvl{14}='batch.keep.rand';tmn{14}=['11:44'];tmf{14}=['15:15'];

pathvl{15}='/oriole5/pu57w52/21209_ACSF/'
catchvl{15}='batch12.keep';tmn{15}=['15:15'];tmf{15}=['21'];

pathvl{16}='/oriole5/pu57w52/21209_ACSF/'
catchvl{16}='batch.keep.rand';tmn{16}=['7'];tmf{16}=['10:52'];

pathvl{17}='/oriole5/pu57w52/21409_APV_5mM/'
catchvl{17}='batch14.keep.rand';tmn{17}=['10:52'];tmf{17}=['14:35'];

pathvl{18}='/oriole5/pu57w52/21409_ACSF/'
catchvl{18}='batch14.keep.rand';tmn{18}=['14:35'];tmf{18}=['21'];

pathvl{19}='/oriole5/pu57w52/21709_ACSF/'
catchvl{19}='batch.keepb';tmn{19}=['12'];tmf{19}=['15:49'];

pathvl{20}='/oriole5/pu57w52/21709_APV_2mM/'
catchvl{20}='batch.keep';tmn{20}=['15:49'];tmf{20}=['19:04'];

pathvl{21}='/oriole5/pu57w52/21709-2/'
catchvl{21}='batch.keep';tmn{21}=['19:04'];tmf{21}=['21'];

pathvl{22}='/oriole5/pu57w52/218wnoff/'
catchvl{22}='batch19.keep.rand';tmn{22}=['7'];tmf{22}=['10:17'];

pathvl{23}='/oriole5/pu57w52/219_APV_2mM/'
catchvl{23}='batch.keep.rand';tmn{23}=['10:17'];tmf{23}=['13:47'];

pathvl{24}='/oriole5/pu57w52/219_ACSF/'
catchvl{24}='batch.keep.rand';tmn{24}=['13:47'];tmf{24}=['21'];

pathvl{25}='/oriole5/pu57w52/21209_ACSF/'
catchvl{25}='batch13';tmn{25}=['7'];tmf{25}=['21'];

pathvl{26}='/oriole5/pu57w52/21409_ACSF/'
catchvl{26}='batch16';tmn{26}=['7'];tmf{26}=['21'];

% pathvl{27}='/oriole5/pu57w52/219_ACSF/'
% catchvl{27}='batch.keep.rand';tmn{24}=['13:47'];tmf{24}=['21'];


usex(1:49)=0;
ac(49)=[0];
avls.analind=[1:49]

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
avls.pitchalign{3}=0;



NFFT(1)=512
NFFT(2)=512
% NFFT(3)=256; 

fbins{1}=[6500,8400];
fbins{2}=[6400,8400];
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
  
acon=[ 1 3 5 6 7 8 10 12 13 15 16 18 19 21 22 24:26] ;  
muon=[2 4 9  11 14 17 20 23]
revanal{1}{1}=[];
revanal{1}{2}=[]
revanal{2}{1}=[]
revanal{2}{2}=[]
muanal{1}{1}=[2 4 9 11 14 17 20 23]


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
avls.mkfv=[25]
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


