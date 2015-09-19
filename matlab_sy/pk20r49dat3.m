%designed to call pitchsyntaxanal, and then inactivanal.m
sumpath='/oriole6/pk20r49/'
matfilename='pathvals1'
prefix='/oriole6/pk20r49/'
sumclearscript;
%change tbinshift to .015 and run again.
%excluding mu ron on 1/25.


wn(1).tmon{1}={'2008-1-31 7' '2008-2-03 7'}
wn(1).tmon{2}={'2008-2-11 7'}
wn(1).tmoff{1}={'2008-2-03 7' '2008-2-7 7'}
wn(1).tmoff{2}={'2008-2-17 7'}
wn(2).tmon=wn(1).tmon;
wn(2).tmoff=wn(1).tmoff;

wnrev(1).tmon{1}={'2008-2-07 7'}
wnrev(1).tmoff{1}={'2008-2-11 7'}
wnrev(2).tmon=wnrev(1).tmon;
wnrev(2).tmoff=wnrev(1).tmoff;


avls.initmean{1}=3245;
avls.initsd{1}=78;

avls.initmean{2}=2178
avls.initsd{2}=137

pathvl{1}='ac250108/'
catchvl{1}='batchcomb'
tmn{1}='7';tmf{1}='13:56'

pathvl{2}='100mu260108real/'
catchvl{2}='batch.keep'
tmn{2}='13:56';tmf{2}='18:51'

pathvl{3}='ac260108/'
catchvl{3}='batchcomb'
tmn{3}='18:51';tmf{3}='21'

pathvl{4}='ac260108/'
catchvl{4}='batch.keep27'
tmn{4}='7';tmf{4}='10:06'

pathvl{5}='200mu270108/'
catchvl{5}='batch.keep'
tmn{5}='10:06:00';tmf{5}='14:17:00'

pathvl{6}='ac270108/'
catchvl{6}='batch27.keep'
tmn{6}='14:17';tmf{6}='21'

pathvl{7}='wnon/'
catchvl{7}='batch29.catch.keep'
tmn{7}='7';tmf{7}='16:32'

pathvl{8}='200muampon12908/'
catchvl{8}='batch29.catch'
tmn{8}='16:32';tmf{8}='19:55'

%there's more data here before switch.
pathvl{9}='ac130/'
catchvl{9}='batch30comb'
tmn{9}='7';tmf{9}='15:06'

pathvl{10}='200muampoff13008/'
catchvl{10}='batch.keep';tmn{10}='15:06:00';tmf{10}='19:03:00'

pathvl{11}='acampoff_300108/'
catchvl{11}='batch.keep'
tmn{11}='19:03';tmf{11}='21';

pathvl{12}='acampoff_300108/'
catchvl{12}='batchcomb'
tmn{12}='7';tmf{12}='15:41';

pathvl{13}='200muampon1.31.08/'
catchvl{13}='batch.catch'
tmn{13}='15:41:00';tmf{13}='19:41:00'

pathvl{14}='acampon13108-3/'
catchvl{14}='batch31'
tmn{14}='19:41';tmf{14}='21';

%this is 2/2 post, where is pre data? BELOW
pathvl{15}='acsf_ampon2/'
catchvl{15}='batch.train.catch'  
tmn{15}='16:34';tmf{15}='21';

pathvl{16}='200muampon202008/'
catchvl{16}='batch.train.catch'
tmn{16}='12:33';tmf{16}='16:34'

pathvl{17}='acampon13108-4/'
catchvl{17}='batch02.catch'
tmn{17}='7';tmf{17}='12:33';

%changed probes 2/4...discount brief muscimol run.

pathvl{18}='probein20408/'
catchvl{18}='batch05.catch' 
tmn{18}='7';tmf{18}='10:48';

pathvl{19}='200mu20508/'
catchvl{19}='batch.catch'
tmn{19}='10:48';tmf{19}='15:05'

pathvl{20}='ac020508/'
catchvl{20}='batch.catch.keep'
tmn{20}='15:05';tmf{20}='21';

pathvl{21}='ac20608temp5/'
catchvl{21}='batch'
tmn{21}='7';tmf{21}='15:12';

pathvl{22}='200mu20608/'
catchvl{22}='batch'
tmn{22}='15:12';tmf{22}='19:28'

pathvl{23}='ac20608-2/'
catchvl{23} ='batch06.catch'
tmn{23}='19:28';tmf{23}='21';

% pathvl{24}='acampon210/'
% catchvl{24}='batch.catch.keep12'
% tmn{24}='7';tmf{24}='10:35';
% 
% pathvl{25}='300mu21208/'
% catchvl{25}='batch.catch'
% tmn{25}='10:35';tmf{25}='14:31'
% 
% pathvl{26}='ac21208ampon/'
% catchvl{26}='batch.catch'
% tmn{26}='14:31';tmf{26}='21'

pathvl{27}='ac21208ampon/'
catchvl{27}='batch13.catch'
tmn{27}='7';  ; tmf{27}='11:25'

pathvl{28}='300mu21308/'
catchvl{28}='batch.catch'
tmn{28}='11:25'; tmf{28}='15:50'

pathvl{29}='ac21308/'
catchvl{29}='batchcomb13'
tmn{29}='15:54';tmf{29}='21';

% pathvl{30}='acsf21508_2/'
% catchvl{30}='batch.catch'
% tmn{30}='7';tmf{30}='13:01';
% 
% pathvl{31}='300mu21608/'
% catchvl{31}='batch.catch'
% tmn{31}='13:36';tmf{31}='19:00'
% 
% pathvl{32}='ac21608/'
% catchvl{32}='batch16.catch'
% tmn{32}='19:00'; tmf{32}='21';

pathvl{33}='ac223/'
catchvl{33}='batch.catch.keep'
tmn{33}='7';tmf{33}='15:54'

pathvl{34}='mu223/'

catchvl{34}='batch'
tmn{34}='15:54'; tmf{34}='18:30';

pathvl{35}='350mu224/'
catchvl{35}='batch.keep'
tmn{35}='7'; tmf{35}='21';

pathvl{36}='200muampon12908/'
catchvl{36}='batch29'
tmn{36}='19:55';tmf{36}='21:00'

pathvl{37}='ampoff_feb8/'
catchvl{37}='batch08.keep'
tmn{37}='7';tmf{37}='13:07'

pathvl{38}='200mu20208/'
catchvl{38}='batch'
tmn{38}='13:07';tmf{38}='17:00'

pathvl{39}='acsf_feb8/'
catchvl{39}='batch08.keep.rand'
tmn{39}='17:00';tmf{39}='21:00'

pathvl{40}='acsf_feb8/'
catchvl{40}='batch09.keep'
tmn{40}='7';tmf{40}='13:20'

pathvl{41}='300mu20908/'
catchvl{41}='batchcomb.keep'
tmn{41}='13:20';tmf{41}='17:35'

pathvl{42}='ac20908/'
catchvl{42}='batch09.keep'
tmn{42}='17:35';tmf{42}='21:00'

pathvl{43}='acampon210/'
catchvl{43}='batch.catch.keep12'
tmn{43}='7';tmf{43}='10:30'

pathvl{44}='300mu21208/'
catchvl{44}='batch.catch'
tmn{44}='10:30';tmf{44}='14:30'

pathvl{45}='ac21208ampon/'
catchvl{45}='batch12.catch.keep.rand'
tmn{45}='14:30';tmf{45}='21:00'

pathvl{46}='fromevren/ac21308/'
catchvl{46}='batch.train.catch'
tmn{46}='7';tmf{46}='11:30'

pathvl{47}='fromevren/mu21408/'
catchvl{47}='batch.train.catch'
tmn{47}='11:30';tmf{47}='15:45'

pathvl{48}='fromevren/ac21408/'
catchvl{48}='batch.train.catch'
tmn{48}='15:45';tmf{48}='21:00'



% pathvl{46}='ac130/'
% catchvl{46}='batch30comb'
% tmn{46}='7';tmf{46}='15:20'
% 
% pathvl{47}='200muampoff13008/'
% catchvl{47}='batch.keep'
% tmn{47}='15:25';tmf{47}='19:02'
% 
% pathvl{48}='acampoff_300108/'
% catchvl{48}='batch.keep.rand'
% tmn{48}='19:02';tmf{48}='21'


pathvl{49}='acsf21508_2/'
catchvl{49}='batch.catch'
tmn{49}='7';tmf{49}='13:25'

pathvl{50}='300mu21608/'
catchvl{50}='batch.catch'
tmn{50}='13:25';tmf{50}='19:00'

pathvl{51}='ac21608/'
catchvl{51}='batch16.catch'
tmn{51}='19:00';tmf{51}='21'


pathvl{52}='ac130/'
catchvl{52}='batch29'
tmn{52}='19:55';tmf{52}='21'



pathvl{53}='ac223-2/'
catchvl{53}='batch23'
tmn{53}='19:50';tmf{53}='21:00'



pathvl{54}='acampon13108-3/'
catchvl{54}='batch01.catch'
tmn{54}='7';tmf{54}='21'


pathvl{55}='probein20408/'
catchvl{55}='batch04.catch'
tmn{55}='7';tmf{55}='21:00'

pathvl{56}='acampon210/'
catchvl{56}='batch11.catch'
tmn{56}='7';tmf{56}='21:00'

pathvl{57}='ac21508_newth/'
catchvl{57}='batch15.catch'
tmn{57}='7';tmf{57}='21:00'

pathvl{58}='ampoff_feb8/'
catchvl{58}='batch07.rand.keep'
tmn{58}='7';tmf{58}='21:00'

pathvl{59}='acsf_ampon3/'
catchvl{59}='batch03.catch.keep'
tmn{59}='7';tmf{59}='21:00'

pathvl{60}='ac20908/'
catchvl{60}='batch10.keep'
tmn{60}='7';tmf{60}='21:00'


% pathvl{36}='probein/'
% catchvl{36}='batchcomb'
% tmn{36}='7'; tmf{36}='21';

numnotes=1
notes='a'

fbins={};
avls.TBINVAR=1;
tbinshft{1}(7:60)=0.024;
tbinshft{1}(1:6)=0.02;
tbinshft{2}(1:60)=0.015
NFFT(1)=256
NFFT(2)=256

fbins{1}=[2700,3700];
fbins{2}=[1700, 2300]
avls.edges{1}=[2700:100:3700];
avls.edges{2}=[1700:100:2300];
clear NT    
NT{1}='z'
NT{2}='a'


PRENT{1}='';PSTNT{1}='';
PRENT{2}='';PSTNT{2}='';
avls.changenote=0;
  %%%plotvals
muon=[2 5 8 10 13 16  19 22  28  34  38 41 44 50 47 ];  
acon=[1 3 4 6 7 9  11 12 14 15 17 18 20 21  23  27 29 33 35  37 39 40 42 43 45 46 48  49 51 52 53 54 55 56 57 58 59 60]
diron=[];

avls.diron=[];
avls.mkfv=[34]
avls.muanal{1}{1}=[2 5 8 10 13 16  19 22  28  34  ];
avls.muanal{1}{2}=[44 50 47]
avls.revanal{1}{1}=[22 38]
avls.maxmeans=[3477];
avls.bname='pk20r49'
graphvals.wn=wn;
graphvals.wnrev=wnrev;
avls.exclude=[0 1]
for ii=1:length(fbins)
    edges{ii}=[fbins{ii}(1):50:fbins{ii}(end)]
end


avls.usex=[0]

clear colvals
colvals{1}=muon;
colvals{2}=acon
colvals{3}=diron;

graphvals.numcol=3
graphvals.col='rkc'
% graphvals.plttext=1


avls.tmon=tmn;
avls.tmoff=tmf;
graphvals.tmn=tmn;
graphvals.tmf=tmf;
graphvals.acon=acon;

avls.edges=edges;

graphvals.colvals=colvals
graphvals.chunkdata=1

if exist('prefix')
    [pathvl]=addprefix(pathvl,prefix);
end
avls.pvls=pathvl
avls.datfile=matfilename;
avls.cvl=catchvl
avls.sumpath=sumpath
avls.mtflnm=matfilename  
avls.NT=NT
avls.acon=acon;
avls.muon=muon;
avls.NFFT=NFFT
avls.fbins=fbins
avls.tshft=tbinshft
avls.numnotes=numnotes
avls.PRENT=PRENT
avls.PSTNT=PSTNT
avls.deadtm=.16/24;
avls.muoffset=1.5/24;
avls.acoffset=1/24;
avls.pretm=5/24

strcmd=['cd ' sumpath 'datasum']
eval(strcmd)

strcmd=['save ' matfilename '.mat avls graphvals'];
eval(strcmd);



   