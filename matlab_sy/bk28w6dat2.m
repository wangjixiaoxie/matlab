%designed to call pitchsyntaxanal, and then inactivanal.m


%excluding mu ron on 1/25.



sumpath='/oriole6/bk28w6/'
matfilename='pathvals1'
prefix='/oriole6/bk28w6/'
sumclearscript;

avls.initmean{1}=4945;
avls.initsd{1}=160;

avls.initmean{2}=6498
avls.initsd{2}=108;


% wn(1).tmon{1}={'2008-7-12 7' '2008-7-16 7'}
% wn(1).tmon{2}={'2008-7-27 7' '2008-7-30 7'}  
% wn(1).tmon{3}={'2008-9-19 7' }
% wn(1).tmoff{1}={'2008-7-16 7' '2008-7-22 7'}
% wn(1).tmoff{2}={'2008-7-30 7' '2008-8-3 7'} 
% wn(1).tmoff{3}={'2008-9-26 7'}
% wn(2).tmon=wn(1).tmon;
% wn(2).tmoff=wn(1).tmoff;
% 
% wn(3).tmon=wn(1).tmon;
% wn(3).tmoff=wn(1).tmoff;

wn(1).tmon{1}= {'2007-11-11 7'}
wn(1).tmoff{1}={'2007-11-24 7'}

wnrev(1).tmon{1}={};
wnrev(2).tmoff{1}={};
graphvals.wnrev=wnrev;
pathvl{1}='screen/'
catchvl{1}='batchcomb'
tmn{1}='7'; tmf{1}='21';

pathvl{2}='temptest2/'
catchvl{2}='batch.keep.rand'
tmn{2}='7';tmf{2}='21';

pathvl{3}='wnon/'
catchvl{3}='batchcomb'
tmn{3}='7';tmf{3}='21';

pathvl{4}='postimplant/'
catchvl{4}='batch.keep.catch'
tmn{4}='7'; tmf{4}='21'

pathvl{5}='mu20010907-2/'
catchvl{5}='batch.keep'
tmn{5}='15:00';tmf{5}='18:58'


% pathvl{6}='probein/'
% catchvl{6}='batch.keep'
% tmn{6}='7';tmf{6}='21';

pathvl{7}='ac20010907/'
catchvl{7}='batch.keep'
tmn{7}='7'; tmf{7}='15:00';

pathvl{8}='ac109007-2/'
catchvl{8}='batch.keep'
tmn{8}='18:58';tmf{8}='21';

 pathvl{9}='ac111007_ampon/'
catchvl{9}='batch13comb'
tmn{9}='7'  ;tmf{9}='9:17'
%on at 817, off at 1355


pathvl{10}='mu111307_ampon/'
catchvl{10}='batchcomb.catch'
tmn{10}='09:17';tmf{10}='13:55'


pathvl{13}='ac111307_ampon/'
catchvl{13}='batch13.keep.catch'
tmn{13}='13:55';tmf{13}='21';

%on at 906, off at 1344


pathvl{14}='mu250_11407_ampon/'
catchvl{14}='batch.catch'
tmn{14}='09:06';tmf{14}='13:44'

% pathvl{15}='/doya2/bk28w6/temp/ac111007_ampon/'
% catchvl{15}='batch13.keep'

pathvl{16}='ac111307_ampon/'
catchvl{16}='batch14.keep'
tmn{16}='7';tmf{16}='09:06'

pathvl{17}='ac111407ampon/'
catchvl{17}='batch2.catch'
tmn{17}='13:44';tmf{17}='21';
%previous_file here is 
pathvl{11}='ac111407ampon/'
catchvl{11}='batch15comb.catch'
tmn{11}='7';tmf{11}='9:30'
% 
% pathvl{18}='mu500_11507_ampon/'
% catchvl{18}='batch.catch'
% tmn{18}='9:30';tmf{18}='13:41'
% 
% pathvl{19}='ac11507ampon/'
% catchvl{19}='batch.keep.catch'
% tmn{19}='13:41';tmf{19}='21'

pathvl{20}='ac11507ampon/'
catchvl{20}='batch16.keep.catch.rand'
tmn{20}='7';tmf{20}='21';

pathvl{21}='ac11507ampon/'
catchvl{21}='batch17.keep.catch.rand'
tmn{21}='7';tmf{21}='21';

pathvl{22}='acprobein1128/'
catchvl{22}='batch.keep.rand'
tmn{22}='7';tmf{22}='21';

pathvl{23}='acprobein1128/'
catchvl{23}='batch30.keep.rand'
tmn{23}='7';tmf{23}='11:04';

pathvl{24}='mu200_113007/'
catchvl{24}='batch.keep.rand'
tmn{24}='11:04';tmf{24}='15:50'

pathvl{25}='ac11307/'
catchvl{25}='batch30'
tmn{25}='15:50';tmf{25}='21';

pathvl{26}='mu500120107/'
catchvl{26}='batch.keep.rand'
tmn{26}='10:55';tmf{26}='14:26'


%pre26
pathvl{12}='ac11307/'
catchvl{12}='batch01'
tmn{12}='7';tmf{12}='10:55'
%these are the extra labels for z and m

%post26
pathvl{27}='ac120107/'
catchvl{27}='batch01comb'
tmn{27}='14:26';tmf{27}='21'


%PRE-ac11007_ampon

% pathvl{28}='mu111307_ampon/'
% catchvl{28}='batch.keep.notcatch'
% tmn{28}='9:17';tmf{28}='14:55';
% 
% %pre
% % pathvl{29}='ac111007_ampon/'
% % catchvl{29}='batch13comb'
% % tmn{29}='7';tmf{29}='9:17';
% 
% %post27
% pathvl{30}='ac111307_ampon/'
% catchvl{30}='batch13.keep.catch'
% tmn{30}='14:55';tmf{30}='21';


%POST ac111307_ampon batch13.keep.catch


%pre acsf= ac1307_ampon %batch14.keep
% pathvl{31}='mu250_11407_ampon/'
% catchvl{31}='batch.keep.notcatch'
% tmn{31}='10:07';tmf{31}='14:44'
% %post acsf=ac111407_ampon batch14comb.catch
% 
% pathvl{32}='ac111307_ampon/'
% catchvl{32}='batch14.keep'
% tmn{32}='7';tmf{32}='10:07'
% 
% pathvl{33}='ac111407ampon/'
% catchvl{33}='batch14comb.catch'
% tmn{33}='14:44';tmf{33}='21'

%preacsf=ac11407_ampon              batch15comb.catch
pathvl{34}='mu500_11507_ampon/'
catchvl{34}='batchcomb'
tmn{34}='9:31';tmf{34}='13:41'

pathvl{35}='ac11507ampon/'
catchvl{35}='batch15.catch'
tmn{35}='13:41';tmf{35}='21'

% pathvl{36}='mu250biotin120607/'
% catchvl{36}='batch06.keep.rand'
% tmn{36}='12:05';tmf{36}='21'
% 
% pathvl{37}='ac120107/'
% catchvl{37}='batch06.keep.rand'
% tmn{37}='7';tmf{37}='12:05'


pathvl{38}='ac111407ampon/'
catchvl{38}='batch15comb.catch'
tmn{38}='7';tmf{38}='10:31'

pathvl{39}='ac111007_ampon/'
catchvl{39}='batch12.catch'
tmn{39}='7';tmf{39}='21'

%postacsf=ac11507ampon batch15.catch

usex(1:29)=0;
avls.analind=[1:29]

numnotes=2
notes='ab'

fbins{1}=[4000,6000];
    
fbins{2}=[6000 8000]

tbinshft{1}=0.092;
tbinshft{2}=.015
NFFT(1)=512
NFFT(2)=512


NT{1}='b'
NT{2}='a'


PRENT{1}='';PSTNT{1}='';
PRENT{2}='-';PSTNT{2}='a'

acon=[1:4  7:9 11 13 16 17  20:23 25 12 27    35  39]
muon=[5 10 14  24 26   34 ]
diron=[];
avls.diron=[];
avls.changenote=[];
avls.mkfv=[39]
avls.muanal{1}{1}=[5 10 14 24 26  34 ];
avls.maxmeans=[3477];
avls.bname='bk28w6'
avls.pretm=4;
assignvals_script;
avls.exclude=[0 0 0]
graphvals.wn=wn;



