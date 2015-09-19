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
matfilename='pathvals3'

extractvals=1

%index into cell array is the note
graphvals.wnon=1;

wnbnds{1}=[3000 4000]
wnbnds{2}=[4000 6000]

% 6/5, 6/6 and 6/7 again 6/14. 6/15, 6/16
% wn{1}.freqlo=[3000 3000 3000 3000 3000     3745 3600 3470 3245]
% wn{1}.freqhi=[3400 3450 3500 3600 3675     4000 4000 4000 4000 ]
% wn{2}.freqlo=[5100 5000  4950 4940 4860     4000 4000 4000 4000]
% wn{2}.freqhi=[6000 6000  6000 6000 6000     4950 4900 4875 5220]
% 
% wn{1}.tmon={'2008-06-05 07:00' '2008-06-05 13:15' '2008-06-05 16:17' '2008-06-06 10:20' '2008-06-07 07:00:00' '2008-06-14 07:00' '2008-06-15 07:00' '2008-06-16 07:00' '2008-06-21 07:00'};
% wn{1}.tmoff=[wn{1}.tmon(2:4) '2008-06-14 07:00' wn{1}.tmon(6:8) '2008-06-16 21:00' '2008-07-02 21:00']
% wn{2}.tmon=wn{1}.tmon;
% wn{2}.tmoff=wn{1}.tmoff;
% 
% wn{3}=wn{2}
% wn{4}=wn{1}
% wn{5}=wn{1};



wn(1).tmon{1}={'2008-6-5 7' '2008-6-8 7'}
wn(1).tmon{2}={'2008-6-19 7' '2008-6-23 7'}  
wn(1).tmoff{1}={'2008-6-8 7' '2008-6-14 7'}
wn(1).tmoff{2}={'2008-6-23 7' '2008-7-08 7'} 
wn(2).tmon=wn(1).tmon;
wn(2).tmoff=wn(1).tmoff;

wn(3).tmon=wn(1).tmon;
wn(3).tmoff=wn(1).tmoff;


wn(4).tmon=wn(1).tmon;
wn(4).tmoff=wn(1).tmoff;

wn(5).tmon=wn(1).tmon;
wn(5).tmoff=wn(1).tmoff;
avls.initmean{1}=3384
avls.initsd{1}=58.5;
avls.maxmeans=[4731 4765]

avls.initmean{2}=5115.5
avls.initsd{2}=95.24


avls.initmean{3}=5115.5
avls.initsd{3}=95.24

avls.initmean{4}=avls.initmean{1}
avls.initsd{4}=avls.initsd{1}

avls.initmean{5}=avls.initmean{4};
avls.initsd{5}=avls.initsd{4};

edges{1}=[6300:55:7500]
edges{2}=[1800:60:2500]
edges{3}=[1800:60:2500]
edges{4}=edges{3}
edges{5}=edges{4}
clear pathvl vals




pathvl{1}='/oriole5/bk15bk14/probein2/'
catchvl{1}='batchcomb'
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
tmn{6}=tmf{5}
tmf{6}='16:29'
pathvl{7}='/oriole5/bk15bk14/500mu6408/'
catchvl{7}='batch'
tmn{7}='16:29'
tmf{7}='19:15'
pathvl{8}='/oriole5/bk15bk14/ac6408-2/'
catchvl{8}='batch06.catch'
tmn{8}='07:00'
tmf{8}='12:45'
pathvl{9}='/oriole5/bk15bk14/ac6408-2/'
catchvl{9}='batch04.keep'
tmn{9}='19:15'
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
tmn{13}='07:00'
tmf{13}='11:11'
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
tmf{22}='14:53'

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
catchvl{30}='batchcomb.keep'
tmn{30}='10:24'
tmf{30}='14:44'

pathvl{31}='/oriole5/bk15bk14/ac622wnon_3/'
catchvl{31}='batch23all.catch.keep'

tmn{31}='07:00';
tmf{31}='13:29';

pathvl{32}='/oriole5/bk15bk14/1000mu623/'
catchvl{32}='batch.catch'
tmn{32}='13:50'
tmf{32}='18:04'
pathvl{33}='/oriole5/bk15bk14/ac623/'
catchvl{33}='batch.keep'
tmn{33}='18:04'
tmf{33}='21:00'

pathvl{34}='/oriole5/bk15bk14/1000mu624/'
catchvl{34}='batchcomb'
tmn{34}='12:16'
tmf{34}='16:16'

pathvl{35}='/oriole5/bk15bk14/ac624/'
catchvl{35}='batch25.keep.catch'
tmn{35}='07:00'
tmf{35}='10:20'

pathvl{36}='/oriole5/bk15bk14/1000mu625/'
catchvl{36}='batch.catch.keep'
tmn{36}='10:20'
tmf{36}='14:24'

pathvl{37}='/oriole5/bk15bk14/ac625/'
catchvl{37}='batch.catch.keep'
tmn{37}='14:24'
tmf{37}='21:00'

pathvl{38}='/oriole5/bk15bk14/ac6808/'
catchvl{38}='batch08.catch.keep'
tmn{38}='15:13'
tmf{38}='21:00'

pathvl{39}='/oriole5/bk15bk14/ac61108/'
catchvl{39}='batch.catch.keep'
tmn{39}='13:30'
tmf{39}='21:00'

pathvl{40}='/oriole5/bk15bk14/ac61308/'
catchvl{40}='batch13.keep.catch'
tmn{40}='16:00'
tmf{40}='21:00'

pathvl{41}='/oriole5/bk15bk14/ac618/'
catchvl{41}='batch18.keep'
tmn{41}='18:37'
tmf{41}='21:00'

pathvl{42}='/oriole5/bk15bk14/ac620/'
catchvl{42}='batch20.keep.rand'
tmn{42}='14:44'
tmf{42}='21:00'

pathvl{43}='/oriole5/bk15bk14/ac624/'
catchvl{43}='batch24.catch.keep'
tmn{43}='16:16'
tmf{43}='21:00'

pathvl{44}='/oriole5/bk15bk14/ac618/'
catchvl{44}='batch20.keep.rand'
tmn{44}='07:00'
tmf{44}='10:24'

pathvl{45}='/oriole5/bk15bk14/ac623/'
catchvl{45}='batch24.catch.keep'
tmn{45}='07:00'
tmf{45}='12:16'


pathvl{46}='/oriole5/bk15bk14/ac61108/'
catchvl{46}='batch13.catch.keep'
tmn{46}='07:00'
tmf{46}='11:44'

pathvl{47}='/oriole5/bk15bk14/ac627/'
catchvl{47}='batch29.keep.catch'
tmn{47}='07:00'
tmf{47}='12:41'

pathvl{48}='/oriole5/bk15bk14/1000mu629/'
catchvl{48}='batch.keep.catch'
tmn{48}='12:41'
tmf{48}='16:44'

pathvl{49}='/oriole5/bk15bk14/ac627/'
catchvl{49}='batch27.catch.keep'
tmn{49}='14:15'
tmf{49}='21:00'

pathvl{50}='/oriole5/bk15bk14/ac625/'
catchvl{50}='batch27.catch.keep'
tmn{50}='07:00'
tmf{50}='10:20'

pathvl{51}='/oriole5/bk15bk14/1000mu627/'
catchvl{51}='batch.catch'
tmn{51}='10:20'
tmf{51}='14:15'

pathvl{52}='/oriole5/bk15bk14/1000mu702/'
catchvl{52}='batch.catch.keep'
tmn{52}='15:10'
tmf{52}='18:39'

pathvl{53}='/oriole5/bk15bk14/ac629/'
catchvl{53}='batch02.catch.keep'
tmn{53}='07:00'
tmf{53}='15:10'

pathvl{54}='/oriole5/bk15bk14/ac702/'
catchvl{54}='batch02.catch'
tmn{54}='18:39'
tmf{54}='21:00'

pathvl{55}='/oriole5/bk15bk14/1000mu704/'
catchvl{55}='batch.catch'
tmn{55}='13:00'
tmf{55}='16:44'

pathvl{56}='/oriole5/bk15bk14/ac702/'
catchvl{56}='batch04.catch.keep'
tmn{56}='7'
tmf{56}='13:00'


pathvl{57}='/oriole5/bk15bk14/ac704/'
catchvl{57}='batch06.catch.keep'
tmn{57}='7'
tmf{57}='12:00'

pathvl{58}='/oriole5/bk15bk14/1000mu706/'
catchvl{58}='batch06.catch.keep'
tmn{58}='12:00'
tmf{58}='16:00'

pathvl{59}='/oriole5/bk15bk14/ac706/'
catchvl{59}='batch06.catch.keep'
tmn{59}='16:00'
tmf{59}='21:00'

pathvl{60}='/oriole5/bk15bk14/ac614/'
catchvl{60}='batch'
tmn{60}='18:38'
tmf{60}='21'



pathvl{61}='/oriole5/bk15bk14/ac629/'
catchvl{61}='batch29.keep.catch'
tmn{61}='16:44'
tmf{61}='21'

pathvl{62}='/oriole5/bk15bk14/ac704/'
catchvl{62}='batch04.keep.catch'
tmn{62}='16:44'
tmf{62}='21'


usex(1:59)=0;
ac(59)=[0];
avls.analind=[1:59]

timon{59}=[];
timoff{59}=[];
date{59}='';
diron=[];
numnotes=2
notes='ab'

fbins={};
tbinshft{1}=0.012;
  tbinshft{2}=0.012;
tbinshft{3}=0.012;
  tbinshft{4}=0.018;
  tbinshft{5}=0.018
NFFT(1)=512
NFFT(2)=512
 NFFT(3)=512
NFFT(4)=512
 NFFT(5)=512

fbins{1}=[3000,4000];
fbins{2}=[4000,6000];
fbins{3}=[4000,6000];
fbins{4}=[3000,4000];
fbins{5}=[3000,4000];

clear NT    
NT{1}='a'
 NT{2}='b'
NT{3}='b'
NT{4}='a'
NT{5}='a'
% NT{3}='g'

PRENT{1}='-';PSTNT{1}='a';
PRENT{2}='';PSTNT{2}='b';
PRENT{3}='b';PSTNT{3}='';
PRENT{4}='a';PSTNT{4}='a';
PRENT{5}='b';PSTNT{5}='a';
  %%%plotvals
acon=[ 1 3 4 6 8 9 11 13 15 16 17 19 22 24 25 27 28 31 33 35 37 38 39 40 41 42 43 44 45 46 47 49 50 53 54 56 57 59 60 61 62];  
muon=[2 5 7 10 12 14 18 20 21 23 26 29 30  32 34 36 20 48 51 52 55 58]
diron=[] 
extraon=[]
avls.muanal{1}{1}=muon;
avls.muanal{2}{1}=muon;
avls.muanal{3}{1}=muon;
avls.muanal{4}{1}=muon;
avls.muanal{5}{1}=muon;

avls.deadtm=.0833/24;
avls.muoffset=1.5/24;
avls.acoffset=1.5/24;
avls.pretm=4/24;
avls.edges=edges;
avls.bname='bk15bk14'

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
avls.mkfv=[32]
avls.exclude=[0 0 0 0 0]
avls.changenote=0;
avls.PRENT=PRENT
avls.PSTNT=PSTNT
avls.repeatanal=[0 0 0 0 0]
avls.bnds{1}='2007-09-20 07:00:00'
avvls.bnds{2}='2007-11-20 07:00:00'
avls.muoffset=1/24;
avls.offset=1;
avls.acoffset=0.5/24;
avls.deadtm=.33/24;
avls.acon=acon;
avls.muon=muon;
avls.diron=diron;

strcmd=['cd ' sumpath 'datasum']
eval(strcmd)

strcmd=['save ' matfilename '.mat avls graphvals'];
eval(strcmd);


