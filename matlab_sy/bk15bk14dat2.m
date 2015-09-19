%designed to call pitchsyntaxanal, and then inactivanal.m
clear date
clear avls
clear ac
clear graphvals
sumpath='/oriole5/bk15bk14/'
strcmd=['cd ' sumpath]
eval(strcmd);
if (~exist('datasum','dir'))
    !mkdir datasum
end
matfilename='pathvals2'

extractvals=1

%index into cell array is the note
wnon=1;
wnbnds{1}=[3000 4000]
wnbnds{2}=[4000 6000]

%6/5, 6/6 and 6/7.
% wn{1}.freq1=[3400 3450 3500 3600 3675     3745 3600 3470]
% wn{2}.freq2=[5100 5000 4950 4900 4875]
% wn{1}.tmon={'2008-06-05 07:00' '2008-06-05 13:15' '2008-06-05 16:17' '2008-06-06 10:20' '2008-06-07 07:00:00'}
% wn{1}.tmoff=[bel_hits{1}.tmon(2:5) '2008-06-07 21:00']
% wn{2}.freq=[5100 5000 4950 4900 4875]
% wn{2}.tmoff=bel_hits{1}.date;
% 
%6/14, 6/15, 6/16
% abov_hits{1}.freq=[3745 3600 3470]
% abov_hits{1}.tmon=['2008-06-14 07:00' '2008-06-15 07:00' '2008-06-16 07:00']
% abov_hits[abov_hits{1}.tmon(2:3) '2008-06-16 21:00']
% bel_hits{2}.freq=[4795 4930 4990]
% bel_hits{2}.date=abov_hits{1}.date

clear pathvl vals




pathvl{1}='/oriole5/bk15bk14/probein2/'
catchvl{1}='batch03.keep.rand'
tmn{1}='07:00'
tmf{1}='11:49'

% timon{1}='13:57:00'
% timoff{1}='18:51:00'
% date{1}='2008-01-26'

pathvl{2}='/oriole5/bk15bk14/500mu6308/'
catchvl{2}='batch.keep'
tmn{2}='11:49'
tmf{2}='15:06'
pathvl{3}='/oriole5/bk15bk14/ac5308/'
catchvl{3}='batch.keep.rand'
tmn{3}='15:06'
tmf{3}='21:00'

pathvl{4}='/oriole5/bk15bk14/ac5308/'
catchvl{4}='batchcomb'
tmn{4}='07:00'
tmf{4}='10:16'
pathvl{5}='/oriole5/bk15bk14/lid6408/'
catchvl{5}='batch'
tmn{5}='10:16'
tmf{5}='13:32'
pathvl{6}='/oriole5/bk15bk14/ac5408/'
catchvl{6}='batch.keep'
tmn{6}=tmn{5}
tmf{6}='21:00'
pathvl{7}='/oriole5/bk15bk14/500mu6408/'
catchvl{7}='batch'
tmn{7}='16:29'
tmf{7}='19:15'
pathvl{8}='/oriole5/bk15bk14/ac6408-2/'
catchvl{8}='batch06.catch'
tmn{8}=tmf{7}
tmf{8}='21:00'
pathvl{9}='/oriole5/bk15bk14/ac6408-2/'
catchvl{9}='batch04.keep'
tmn{9}='7:00'
tmf{9}='21:00'

pathvl{10}='/oriole5/bk15bk14/500mu6608/'
catchvl{10}='batch.catch'
tmn{10}='12:45'
tmf{10}='17:05'
pathvl{11}='/oriole5/bk15bk14/ac5608/'
catchvl{11}='batch.catch'
tmn{11}=tmf{10}
tmf{11}='21:00'
pathvl{12}='/oriole5/bk15bk14/500mu6808/'
catchvl{12}='batch.catch.keep'
tmn{12}='11:11'
tmf{12}='15:13'
pathvl{13}='/oriole5/bk15bk14/ac5708newtemp/'
catchvl{13}='batch08.catch.keep'
tmn{13}='15:13'
tmf{13}='21:00'
pathvl{14}='/oriole5/bk15bk14/1000mu6908/'
catchvl{14}='batch.catch'
tmn{14}='11:21'
tmf{14}='16:01'
pathvl{15}='/oriole5/bk15bk14/ac6808/'
catchvl{15}='batch09.catch.keep.rand'
tmn{15}='07:00'
tmf{15}='11:21'
pathvl{16}='/oriole5/bk15bk14/ac6908/'
catchvl{16}='batch.catch.keep'
tmn{16}='16:01'
tmf{16}='21:00'

pathvl{17}='/oriole5/bk15bk14/ac6908/'
catchvl{17}='batch11.catch.keep'
tmn{17}='07:00'
tmf{17}='9:39'

pathvl{18}='/oriole5/bk15bk14/1000mu61108/'
catchvl{18}='batch.catch.keep'
tmn{18}='09:39'
tmf{18}='13:30'
pathvl{19}='/oriole5/bk15bk14/ac61408contingswitch/'
catchvl{19}='batch.catch'
tmn{19}='07:00'
tmf{19}='15:17'

pathvl{20}='/oriole5/bk15bk14/1000mu61308/'
catchvl{20}='batch.catch.keep'
tmn{20}='11:44'
tmf{20}='16:00'
pathvl{21}='/oriole5/bk15bk14/1000mu61408/'
catchvl{21}='batch.catch.keep'
tmn{21}='15:17'
tmf{21}='18:38'

pathvl{22}='/oriole5/bk15bk14/ac614_newth/'
catchvl{22}='batch.catch.keep'
tmn{22}='07:00'
tmf{22}='18:11'

pathvl{23}='/oriole5/bk15bk14/1000mu615/'
catchvl{23}='batch.catch.keep'
tmn{23}='14:53'
tmf{23}='18:11'

pathvl{24}='/oriole5/bk15bk14/ac615/'
catchvl{24}='batch.catch.keep'
tmn{24}='18:11'
tmf{24}='21:00'

pathvl{25}='/oriole5/bk15bk14/ac615/'
catchvl{25}='batch17.keep.rand'
tmn{25}='07:00'
tmf{25}='11:59'

pathvl{26}='/oriole5/bk15bk14/1000mu617/'
catchvl{26}='batch.keep'
tmn{26}='11:59'
tmf{26}='15:53'

pathvl{27}='/oriole5/bk15bk14/ac617/'
catchvl{27}='batch17.keep.rand'
tmn{27}='15:53'
tmf{27}='21:00'

pathvl{28}='/oriole5/bk15bk14/ac617/'
catchvl{28}='batch18.keep'
tmn{28}='7:00'
tmf{28}='14:27'

pathvl{29}='/oriole5/bk15bk14/1000mu618/'
catchvl{29}='batch.keep'
tmn{29}='14:27'
tmf{29}='18:37'

pathvl{30}='/oriole5/bk15bk14/1000mu620/'
catchvl{30}='batch'
tmn{30}='10:24'
tmf{30}='14:44'

pathvl{31}='/oriole5/bk15bk14/ac622wnon_3/'
catchvl{31}='batch23.keep.catch'

tmn{31}='07:00';
tmf{31}='13:29';

pathvl{32}='/oriole5/bk15bk14/1000mu623/'
catchvl{32}='batch.catch'

pathvl{33}='/oriole5/bk15bk14/ac623/'
catchvl{33}='batch.keep'


pathvl{34}='/oriole5/bk15bk14/1000mu624/'
catchvl{34}='batch.catch'

pathvl{35}='/oriole5/bk15bk14/1000mu625/'
catchvl{35}='batch.catch'


tmn{35}=[];
tmf{35}=[];
timon{35}=[];
timoff{35}=[];
date{35}='';

usex(1:35)=0;
ac(35)=[0];
avls.analind=[1:35]

numnotes=2
notes='ab'

fbins={};
tbinshft{1}=0.016;
  tbinshft{2}=0.015;

NFFT(1)=512
NFFT(2)=1024
 

fbins{1}=[3000,4000];
fbins{2}=[4000,6000];

clear NT    
NT{1}='a'
 NT{2}='b'

% NT{3}='g'

PRENT{1}='--';PSTNT{1}='';
PRENT{2}='';PSTNT{2}='';

  %%%plotvals
acon=[ 1 3 4 6 8 9 11 13 15 16 17 19 22 24 25 27 28 31 33 ];  
muon=[2 5 7 10 12 14 18 20 21 23 26 29 30  32 34 35]
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

% graphvals.wnon=wnon;
% graphvals.ahits=above_hits;
% graphvals.bhits=bel_hits


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
avls.mkfv=[34:35]
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


