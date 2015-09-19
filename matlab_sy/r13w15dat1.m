%designed to call pitchsyntaxanal, and then inactivanal.m
clear date
clear avls
clear ac
clear graphvals
sumpath='/oriole1/r13w15/'
strcmd=['cd ' sumpath]
eval(strcmd);
if (~exist('datasum','dir'))
    !mkdir datasum
end
matfilename='pathvals1'

extractvals=1


% % white noise information
% wnon{1}=  '2008-04-29  07:00:00'
% minbnds{1}=6900
% maxbnds{1}=7430
% wnoff{1} ='2008-5-14 07:00:00'
% 
% wnon{2}='2008-03-12 07:00:00'
% wnoff{2}='2008-03-14 07:00:00'
% minbnds{2}=6700
% maxbnds{2}=7450
% 
% wnon{3}='2008-03-14 07:00:00'
% wnoff{3}='2008-03-17 07:00:00'
% minbnds{3}=6700
% maxbnds{3}=7630
% 
% 
% wnon{4}='2008-03-19 14:00:00'
% wnoff{4}='2008-03-22 07:00:00'
% minbnds{4}=7280 
% maxbnds{4}=8000
% 
% wnon{5}='2008-03-22 7:00:00'
% wnoff{5}='2008-03-26 07:00:00'
% minbnds{5}=6980
% maxbnds{5}=8000




clear pathvl vals




pathvl{1}='/oriole1/r13w15/probout3/'
catchvl{1}='batch17'
tmn{1}='07:00'
tmf{1}='11:49'

% timon{1}='13:57:00'
% timoff{1}='18:51:00'
% date{1}='2008-01-26'

pathvl{2}='/oriole1/r13w15/lid617/'
catchvl{2}='batch.keep.rand'
tmn{2}='11:49'
tmf{2}='15:06'
pathvl{3}='/oriole1/r13w15/ac617/'
catchvl{3}='batch.keep'
tmn{3}='15:06'
tmf{3}='21:00'

pathvl{4}='/oriole1/r13w15/ac618/'
catchvl{4}='batch19.keep'
tmn{4}='07:00'
tmf{4}='10:16'
pathvl{5}='/oriole1/r13w15/lid619/'
catchvl{5}='batch.keep.rand'
tmn{5}='10:16'
tmf{5}='13:32'
pathvl{6}='/oriole1/r13w15/ac619/'
catchvl{6}='batch'
tmn{6}=tmn{5}
tmf{6}='21:00'


pathvl{7}='/oriole1/r13w15/300mu620/'
catchvl{7}='batch'
tmn{6}=tmn{5}
tmf{6}='21:00'

pathvl{8}='/oriole1/r13w15/ac620/'
catchvl{8}='batchcomb'

pathvl{9}='/oriole1/r13w15/ac619/'
catchvl{9}='batch20.keep.rand'

pathvl{10}='/oriole1/r13w15/ac622wnon/'
catchvl{10}='batch23.catch.keep'

pathvl{10}='/oriole1/r13w15/ac622wnon/'
catchvl{10}='batch23.catch.keep'
timon{10}='';
timoff{10}='';
date{10}='';

pathvl{11}='/oriole1/r13w15/300mu623/'
catchvl{11}='batch'

pathvl{12}='/oriole1/r13w15/ac623/'
catchvl{12}='batch.catch.keep'

pathvl{13}='/oriole1/r13w15/500mu624/'
catchvl{13}='batch'

pathvl{14}='/oriole1/r13w15/ac624/'
catchvl{14}='batch.keep.rand'

pathvl{15}='/oriole1/r13w15/1000biot625/'
catchvl{15}='batch.keep'

timon{15}='';
timoff{15}='';
date{15}='';


usex(1:15)=0;
ac(15)=[0];
avls.analind=[1:15]

numnotes=1
notes='a'

fbins={};
tbinshft{1}=0.04;
 

NFFT(1)=1024
 

fbins{1}=[6000,8000];


clear NT    
NT{1}='a'
 

% NT{3}='g'

PRENT{1}='--';PSTNT{1}='';
PRENT{2}='';PSTNT{2}='';

  %%%plotvals
acon=[ 1 3 4 6 8 9 10 12 14  ];  
muon=[2 5 7 11 13 15]
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

graphvals.timon=timon;
graphvals.timoff=timoff;
graphvals.acon=ac;
graphvals.date=date;

graphvals.wnon=wnon;
graphvals.wnoff=wnoff;
graphvals.minbnds=minbnds;
graphvals.maxbnds=maxbnds;


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
avls.mkfv=[15]
avls.PRENT=PRENT
avls.PSTNT=PSTNT
avls.repeatanal=[0 0 0]
avls.bnds{1}='2007-09-20 07:00:00'
avvls.bnds{2}='2007-11-20 07:00:00'
avls.muoffset=1.25/24;
avls.offset=1;
avls.acoffset=0.75/24;
avls.deadtm=.2/24;avls
avls.acon=acon;
avls.muon=muon;

strcmd=['cd ' sumpath 'datasum']
eval(strcmd)

strcmd=['save ' matfilename '.mat avls graphvals'];
eval(strcmd);


