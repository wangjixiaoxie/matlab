%master script, load and combine various fvs, from postscreen./50microstim,
%and 10 microstim

% pathvl{1}='/doyale/twarren/r95pk42/templtest/'
% catchvl{1}='batch.catch.rand'
% pathvl{2}='/doyale/twarren/r95pk42/r95pk42screen/'
% catchvl{2}='batch2.keep.rand'
% pathvl{3}='/doyale/twarren/r95pk42/wnon/'
% catchvl{3}='batch.keep.catch.rand'

%need to organize and label wnon3, where the stable run is located


% pathvl{1}='/doyale/twarren/r95pk42/templtest/'
% catchvl{1}='batch.catch.rand'
% pathvl{2}='/doyale/twarren/r95pk42/r95pk42screen/'
% catchvl{2}='batch2.keep.rand'
clear avls
avls.baspath='/oriole7/dir1/r95pk42'

pvl{1}='/templteste2';
cvl{1}='batch1920comb'
pvl{2}='/templteste3'
cvl{2}='batch.catch.keep'
pvl{3}='/wnon2'
cvl{3}='batch2223'
pvl{4}='/wnon2'
cvl{4}='batch24.catch.keep'
pvl{5}=pvl{4}
cvl{5}='batch25.catch.keep'
pvl{6}='/26files'
cvl{6}='batch26.keep.catch'
pvl{7}='/27files'
cvl{7}='batch27.keep.catch.rand'
pvl{8}='/30files'
cvl{8}='batch30.keep.catch'
pvl{9}='/wnon3'
cvl{9}='batch31.keep.catch'
pvl{10}='/wnon3/01files'
cvl{10}='batch01.keep.catch'
pvl{11}=pvl{9}
cvl{11}='batch02.keep.catch'

pvl{12}='/templatest';
cvl{12}='batch.keep.rand'
pvl{13}='/templatest2/'
cvl{13}='batch.keep.rand'
pvl{14}='/templatest3/'
cvl{14}='batch.keep.catch'
pvl{15}='/wnona'
cvl{15}='batch28.keep.catch'
pvl{16}='/wnona/'
cvl{16}='batch30.keep.catch'
pvl{17}=pvl{16}
cvl{17}='batch01.keep.catch'
pvl{18}=pvl{16}
cvl{18}='batch03.keep.catch'
pvl{19}='/wnoffa/'
cvl{19}='batchcomb.rand'

pvl{20}='/synshift2/10files/'
cvl{20}='batch10.keep.catch'
pvl{21}='/synshift2/'
cvl{21}='batch11.keep.catch'
pvl{22}=pvl{21}
cvl{22}='batch12comb'

pvl{23}=pvl{21}
cvl{23}='batch1314a.keep.catch'

pvl{24}=pvl{21}
cvl{24}='batch14b1516a.keep.catch'

pvl{25}='/synshift2wnoff/'
cvl{25}='batchcomb'

pvl{26}='/wnon3-2/06files'
cvl{26}='batch.rand.keep.rand'

pvl{27}='/wnon3-2/08files'
cvl{27}='batch.rand.keep.rand'

pvl{28}='/wnon3-2/10files'
cvl{28}='batch.rand.keep.rand'

pvl{29}='/wnon3-2/09files'
cvl{29}='batch.rand.keep'

pvl{30}='/wnon3-2/07files'
cvl{30}='batch.rand.keep.rand'


%ADD THIRD SYNTAX SHIFT - synshift2, synshift2wnoff??? 


avls.pvls=pvl;
avls.cvl=cvl;
avls.ANALPTNT=1;
avls.SEQNT=1;
avls.DOCONTSYNANAL=0;
avls.DOCONTPTANAL=0;
for ii=1:length(pvl)
avls.EV4(ii)=0;
end
avls.fbins{1}=[3000 3500];
avls.fbins{2}=[2000 3000];
avls.fbins{3}=[2000 3000];

avls.NT{1}='e';
avls.NT{2}='a';
avls.NT{3}='i'

for ii=1:3
    avls.PRENT{ii}=''
end
avls.tbinshft(1)=.025;
avls.tbinshft(2)=.025
avls.tbinshft(3)=.025;



avls.NFFT(1)=256;
avls.NFFT(2)=256;
avls.NFFT(3)=256;
avls.MKFV=[1:30]

avls.PTIND=[20]
avls.SEQIND=[1:30]


avls.bastms{1}={'2007-5-19' '2007-5-21'}
avls.wntms{1}={'2007-5-22' '2007-5-30'},
avls.rectms{1}={'2007-5-31' '2007-6-04'}
avls.offset_dys(1)=[4];
avls.bastms{2}={'2007-6-22' '2007-6-26'}
avls.wntms{2}={'2007-6-27' '2007-7-04'}
avls.rectms{2}={'2007-7-04' '2007-7-06'}
avls.offset_dys(2)=[4];
avls.bastms{3}={'2007-9-10' '2007-9-10'}
avls.wntms{3}={'2007-9-11' '2007-9-16'}
avls.rectms{3}={'2007-9-17' '2007-9-18'}
avls.offset_dys(3)=[4];

avls.stabwin={'2007-06-06' '2007-06-09'}


avls.SEQTRGNT{1}=1;
avls.SEQTRGNT{2}=2;
avls.SEQTRGNT{3}=1;
avls.ALLSEQNT=[1 2 3]
avls.LEARNANALIND=[1 2]

