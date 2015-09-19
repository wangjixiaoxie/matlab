%designed to call pitchsyntaxanal, and then inactivanal.m
clear date
clear avls
sumpath='/oriole/pk32bk28/'
strcmd=['cd ' sumpath]
eval(strcmd);
if (~exist('datasum','dir'))
    !mkdir datasum
end
matfilename='pathvals1'

extractvals=1


%white noise information
wnon{1}=  '2008-03-09  07:00:00'
minbnds{1}=7320
maxbnds{1}=8000
wnoff{1} ='2008-3-12 07:00:00'

wnon{2}='2008-03-12 07:00:00'
wnoff{2}='2008-03-14 07:00:00'
minbnds{2}=7450
maxbnds{2}=8000

wnon{3}='2008-03-14 07:00:00'
wnoff{3}='2008-03-17 07:00:00'
minbnds{3}=7630
maxbnds{3}=8000


wnon{4}='2008-03-19 14:00:00'
wnoff{4}='2008-03-22 07:00:00'
minbnds{4}=6800
maxbnds{4}=7280

wnon{5}='2008-03-22 7:00:00'
wnoff{5}='2008-03-26 07:00:00'
minbnds{5}=6800
maxbnds{5}=6980




clear pathvl vals




pathvl{1}='/oriole/pk32bk28/ac3608b/'
catchvl{1}='batch.keep.rand'
timon=[]
timoff=[]
date{1}=''
% timon{1}='13:57:00'
% timoff{1}='18:51:00'
% date{1}='2008-01-26'

pathvl{2}='/oriole/pk32bk28/300muampon310/'
catchvl{2}='batch'
timon{2}='13:45:00'
timoff{2}='15:31:00'
date{2}='2008-03-10'

pathvl{3}='/oriole/pk32bk28/ac31008ampon/'
catchvl{3}='batch.catch'
acon(3)=2;

pathvl{4}='/oriole/pk32bk28/300mu3808/'
catchvl{4}='batch.keep.rand'
timon{4}='08:49:00'
timoff{4}='12:58:00'
date{4}='2008-03-08'

pathvl{5}='/oriole/pk32bk28/ac3808ampon/'
catchvl{5}='batch.catch.keep'
acon(5)=4;

pathvl{6}='/oriole/pk32bk28/ampon3908probein/'
catchvl{6}='batch.catch'

pathvl{7}='/oriole/pk32bk28/temptest/'
catchvl{7}='batch02.keep.rand'

pathvl{8}='/oriole/pk32bk28/200mu31108/'
catchvl{8}='batch.catch'
timon{8}='13:48:00'
timoff{8}='15:31:00'
date{8}='2008-03-11'

pathvl{9}='/oriole/pk32bk28/ac31108/'
catchvl{9}='batch.catch'
ac(9)=8;

pathvl{10}='/oriole/pk32bk28/ac31208newtemp/'
catchvl{10}='batch.catch'

pathvl{11}='/oriole/pk32bk28/200mu31308/'
catchvl{11}='batch.catch'
timon{11}='11:51:00'
timoff{11}='15:55:00'
date{11}='2008-03-13'

pathvl{12}='/oriole/pk32bk28/ac31308newtemp/'
catchvl{12}='batch.catch'
ac(12)=11;

pathvl{13}='/oriole/pk32bk28/lid314/'
catchvl{13}='batch.catch'
timon{13}='17:11:00'
timoff{13}='19:20:00'
date{13}='2008-03-14'


pathvl{14}='/oriole/pk32bk28/ac31408/'
catchvl{14}='batch.catch.keep'
ac(14)=13;

pathvl{15}='/oriole/pk32bk28/lid315/'
catchvl{15}='batch.catch'
timon{15}='13:28:00'
timoff{15}='15:38:00'
date{15}='2008-03-15'

pathvl{16}='/oriole/pk32bk28/ac315/'
catchvl{16}='batch.catch'
ac(15)=16;

pathvl{17}='/oriole/pk32bk28/ac315newtemp/'
catchvl{17}='batch.catch.keep'

pathvl{18}='/oriole/pk32bk28/ac31608_wnoff/'
catchvl{18}='batch.keep.rand'

pathvl{19}='/oriole/pk32bk28/lid317wnoff/'
catchvl{19}='batch.keep'
timon{19}='08:25:00'
timoff{19}='11:52:00'
date{19}='2008-03-17'

pathvl{20}='/oriole/pk32bk28/ac317/'
catchvl{20}='batch.keep.rand'
ac(20)=19;

pathvl{21}='/oriole/pk32bk28/lid318/'
catchvl{21}='batch.keep'

pathvl{22}='/oriole/pk32bk28/ac318/'
catchvl{22}='batch.keep'
timon{22}='15:07:00'
timoff{22}='17:19:00'
date{22}='2008-01-26'
ac(22)=21;


pathvl{23}='/oriole/pk32bk28/ac318/'
catchvl{23}='batch2.keep.rand'

pathvl{24}='/oriole/pk32bk28/lid319/'
catchvl{24}='batch.keep.rand'
timon{24}='08:25:00'
timoff{24}='11:22:00'
date{24}='2008-03-19'

pathvl{25}='/oriole/pk32bk28/lid319/'
catchvl{25}='batch2.keep'
timon{25}='08:25:00'
timoff{25}='11:22:00'
date{25}='2008-03-19'

pathvl{26}='/oriole/pk32bk28/ac319/'
catchvl{26}='batch.keep.rand'
ac(26)=25;
pathvl{27}='/oriole/pk32bk28/ac319wnon/'
catchvl{27}='batch.catch'

pathvl{28}='/oriole/pk32bk28/lid320ampon/'
catchvl{28}='batch.catch'
timon{28}='16:13:00'
timoff{28}='18:55:00'
date{28}='2008-03-20'

pathvl{29}='/oriole/pk32bk28/ac320ampon/'
catchvl{29}='batch.catch.keep'
ac(29)=28;
pathvl{30}='/oriole/pk32bk28/lid321ampon/'
catchvl{30}='batch.catch'

timon{30}='10:18:00'
timoff{30}='13:37:00'
date{30}='2008-03-21'

pathvl{31}='/oriole/pk32bk28/ac321ampon/'
catchvl{31}='batch'
ac(31)=30;
pathvl{32}='/oriole/pk32bk28/ac321newtemp/'
catchvl{32}='batch23.catch.keep'

pathvl{33}='/oriole/pk32bk28/lid324/'
catchvl{33}='batch.catch.keep'
timon{33}='10:24:00'
timoff{33}='13:59:00'
date{33}='2008-03-24'

pathvl{34}='/oriole/pk32bk28/ac324/'
catchvl{34}='batch.catch.keep'
ac(34)=33;
pathvl{35}='/oriole/pk32bk28/ac323ampon/'
catchvl{35}='batch.catch.keep'

timon{35}=''
timoff{35}=''
date{35}=''
usex(1:35)=0;
ac(35)=0;
avls.analind=[1:35]

numnotes=1
notes='a'

fbins={};
tbinshft{1}=0.035;
  tbinshft{2}=0.02
 tbinshft{3}=0.03
NFFT(1)=1024
 NFFT(2)=512
NFFT(3)=512

fbins{1}=[6000,8000];
fbins{2}=[1900,2500];
 fbins{3}=[1100,1800];

clear NT    
NT{1}='a'
 NT{2}='b'
 NT{3}='c'
% NT{3}='g'

PRENT{1}='';PSTNT{1}='';
PRENT{2}='';PSTNT{2}='';
 PRENT{3}='';PSTNT{3}='';
   
  %%%plotvals
muon=[2 4 8 11 13 15 19 21 24 25 28 30 33];  
acon=[1 3 5 6 7 9 10 12 14 16 17 18 20 22 23 26 27 29 31 32 34 35]
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
graphvals.edges{1}=[6000:60:8000]
graphvals.edges{2}=[1900:30:2500]
graphvals.edges{3}=[1100:30:1800]

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
avls.numnotes=numnotes
avls.mkfv=[]
avls.PRENT=PRENT
avls.PSTNT=PSTNT
avls.repeatanal=[0 0 0]
avls.bnds{1}='2007-09-20 07:00:00'
avvls.bnds{2}='2007-11-20 07:00:00'
avls.offset=1.75;
avls.acoffset=1;

strcmd=['cd ' sumpath 'datasum']
eval(strcmd)

strcmd=['save ' matfilename '.mat avls graphvals'];
eval(strcmd);


