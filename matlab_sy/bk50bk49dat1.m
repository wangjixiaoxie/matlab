%designed to call pitchsyntaxanal, and then inactivanal.m
clear date
clear avls
clear ac
clear graphvals
clear tmn
clear tmf
sumpath='/oriole3/bk50bk49/'
strcmd=['cd ' sumpath]
eval(strcmd);
if (~exist('datasum','dir'))
    !mkdir datasum
end
matfilename='pathvals2'

extractvals=1

%index into cell array is the note
graphvals.wnon=1;

wnbnds{1}=[3000 4000]


% 6/5, 6/6 and 6/7 again 6/14. 6/15, 6/16
wn{1}.freqlo=[3000 ]
wn{1}.freqhi=[3510]

wn{1}.tmon={'2008-06-30 07:00'};
wn{1}.tmoff=['2008-07-04 07:00']

wn{2}=wn{1};
wn{3}=wn{1};

clear pathvl vals




pathvl{1}='/oriole3/bk50bk49/ac628/'
catchvl{1}='batch29.keep.rand'
tmn{1}='07:00'
tmf{1}='12:35'


% timon{1}='13:57:00'2
% timoff{1}='18:51:00'
% date{1}='2008-01-26'

pathvl{2}='/oriole3/bk50bk49/600mu629/'
catchvl{2}='batch.keep.rand'
tmn{2}='12:35'
tmf{2}='16:11'

pathvl{3}='/oriole3/bk50bk49/ac629/'
catchvl{3}='batch.keep.rand'
tmn{3}='16:11';tmf{3}='21:00'

pathvl{4}='/oriole3/bk50bk49/probein/'
catchvl{4}='batch27.keep.rand'
tmn{4}='7';tmf{4}='10:33'

pathvl{5}='/oriole3/bk50bk49/500mu627/'
catchvl{5}='batch.keep.rand'
tmn{5}='10:33';tmf{5}='14:21'

pathvl{6}='/oriole3/bk50bk49/ac627/'
catchvl{6}='batch27.keep.rand'
tmn{6}='14:21';tmf{6}='21'

pathvl{7}='/oriole3/bk50bk49/ac627/'
catchvl{7}='batch28.keep.rand'
tmn{7}='7';tmf{7}='12:55'

pathvl{8}='/oriole3/bk50bk49/750mu628/'
catchvl{8}='batch.keep.rand'
tmn{8}=tmf{7};tmf{8}='15:57'

pathvl{9}='/oriole3/bk50bk49/ac628/'
catchvl{9}='batch28.keep.rand'
tmn{9}=tmf{8};tmf{9}='21'

pathvl{10}='/oriole3/bk50bk49/ac628/'
catchvl{10}='batch29.keep.rand'
tmn{10}='7';tmf{10}='12:35'

pathvl{11}='/oriole3/bk50bk49/wnon630/'
catchvl{11}='batch02.keep.catch'
tmn{11}='7';tmf{11}='10:58'

pathvl{12}='/oriole3/bk50bk49/750mu702/'
catchvl{12}='batchcomb'
tmn{12}='10:58';tmf{12}='15:05'

pathvl{13}='/oriole3/bk50bk49/ac702/'
catchvl{13}='batch02.catch.keep'
tmn{13}='15:05';tmf{13}='21'

pathvl{14}='/oriole3/bk50bk49/ac702/'
catchvl{14}='batch03.catch.keep.rand'
tmn{14}='7';tmf{14}='11:38'

pathvl{15}='/oriole3/bk50bk49/750mu703/'
catchvl{15}='batch.catch.keep'
tmn{15}='11:38';tmf{15}='15:35'

pathvl{16}='/oriole3/bk50bk49/ac703/'
catchvl{16}='batch.catch'
tmn{16}='15:36';tmf{16}='21:00'

pathvl{17}='/oriole3/bk50bk49/ac704/'
catchvl{17}='batch06.catch.keep'
tmn{17}='7';tmf{17}='12:00'

pathvl{18}='/oriole3/bk50bk49/750mu706/'
catchvl{18}='batch06.catch.keep'
tmn{18}='12:00';tmf{18}='15:50'

pathvl{19}='/oriole3/bk50bk49/ac706/'
catchvl{19}='batch06.catch.keep'
tmn{19}='15:50';tmf{19}='21:00'




usex(1:19)=0;
avls.analind=[1:19]

date{19}='';

numnotes=2
notes='kaa'

fbins={};
tbinshft{1}=0.01;
  tbinshft{2}=0.01;
    tbinshft{3}=0.01;
NFFT(1)=512
NFFT(2)=512
NFFT(3)=512 

fbins{1}=[3000,4000];
fbins{2}=[6000,8000];
fbins{3}=[6000,8000];


clear NT    
NT{1}='k'
 NT{2}='a'
NT{3}='a'

PRENT{1}='';PSTNT{1}='';
PRENT{2}='';PSTNT{2}='';
PRENT{3}='k';PSTNT{3}='a';

  %%%plotvals
acon=[ 1 3 4 6 7 9 10 11 13 14 16 17 19];  
muon=[2 5 8 12 15 18];
diron=[] 
extraon=[]

clear colvals
colvals{1}=muon;
colvals{2}=acon


graphvals.numcol=2
graphvals.col='rk'
graphvals.plttext=1
graphvals.txtht=[3260 3200]
graphvals.tickht=8000
graphvals.edges{1}=[6700:80:7900]
graphvals.edges{2}=[2000:40:2500]
graphvals.edges{3}=[1300:40:1800]
graphvals.colvals=colvals
graphvals.pltpnts=1

graphvals.date=date;

graphvals.wn=wn;


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
avls.tmon=tmn;
avls.tmoff=tmf;
avls.date=date;
avls.numnotes=numnotes
avls.mkfv=[17:19]
avls.PRENT=PRENT
avls.PSTNT=PSTNT
avls.repeatanal=[0 0 0 0 0]
avls.bnds{1}='2007-09-20 07:00:00'
avvls.bnds{2}='2007-11-20 07:00:00'
avls.muoffset=1.5/24;
avls.offset=1;
avls.acoffset=1/24;
avls.deadtm=.33/24;avls
avls.acon=acon;
avls.muon=muon;

strcmd=['cd ' sumpath 'datasum']
eval(strcmd)

strcmd=['save ' matfilename '.mat avls graphvals'];
eval(strcmd);


