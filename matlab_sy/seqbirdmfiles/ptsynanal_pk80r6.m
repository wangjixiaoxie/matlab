%master script, load and combine various fvs, from postscreen./50microstim,
%and 10 microstim

% pathvl{1}='/doyale/twarren/r95pk42/templtest/'
% catchvl{1}='batch.catch.rand'
% pathvl{2}='/doyale/twarren/r95pk42/r95pk42screen/'
% catchvl{2}='batch2.keep.rand'
% pathvl{3}='/doyale/twarren/r95pk42/wnon/'
% catchvl{3}='batch.keep.catch.rand'




% pathvl{1}='/doyale/twarren/r95pk42/templtest/'
% catchvl{1}='batch.catch.rand'
% pathvl{2}='/doyale/twarren/r95pk42/r95pk42screen/'
% catchvl{2}='batch2.keep.rand'
clear avls
avls.baspath='/oriole7/dir1/pk80r6'


%START OFF BY LABELING KEY DATA

%pre data=11/20
pvl{1}='/wnon2pre'
cvl{1}='batchcomb'

%seqshift - synshift2  %11-22 to 11-27
pvl{2}='/synshift2/27files'
cvl{2}='batch27.catch.keep'



%recovery data  11-30....
pvl{3}='/wnoff2/30files'
cvl{3}='batch30.catch.keep'


% pvl{3}='/wnon2'
% cvl{3}='batch2223'
% pvl{4}='/wnon2'
% cvl{4}='batch24.catch.keep'
% pvl{5}=pvl{4}
% cvl{5}='batch25.catch.keep'
% pvl{6}='/wnon3'
% cvl{6}='batch26.keep.catch'
% pvl{7}='/wnon3'
% cvl{7}='batch27.keep.catch.rand'
% pvl{8}=pvl{7}
% cvl{8}='batch30.keep.catch'
% pvl{9}=pvl{8}
% cvl{9}='batch31.keep.catch'
% pvl{10}=pvl{9}
% cvl{10}='batch01.keep.catch'
% pvl{11}=pvl{10}
% cvl{11}='batch02.keep.rand'
% 
% pvl{12}='/templatest';
% cvl{12}='batch.keep.rand'
% pvl{13}='/templatest2/'
% cvl{13}='batch.keep.rand'
% pvl{14}='/templatest3/'
% cvl{14}='batch.keep.catch'
% pvl{15}='/wnona'
% cvl{15}='batch28.keep.catch'
% pvl{16}='/wnona/'
% cvl{16}='batch30.keep.catch'
% pvl{17}=pvl{16}
% cvl{17}='batch01.keep.catch'
% pvl{18}=pvl{16}
% cvl{18}='batch03.keep.catch'
% pvl{19}='/wnoffa/'
% cvl{19}='batchcomb.rand'
% 
% % pvl{20}='/synshift2/'
% cvl{20}='batch10.keep.catch'
% pvl{21}=pvl{20}
% cvl{21}='batch11.keep.catch'
% pvl{22}=pvl{20}
% cvl{22}='batch12comb'
% 
% pvl{23}=pvl{20}
% cvl{23}='batch1314a.keep.catch'
% 
% pvl{24}=pvl{20}
% cvl{24}='batch14b1516a.keep.catch'
% 
% pvl{25}='/synshift2wnoff/'
% cvl{25}='batchcomb'
% 
% 


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
avls.fbins{4}=[2000 3000];

avls.NT{1}='c';
avls.NT{2}='e';
avls.NT{3}='z'
avls.NT{4}='i'

for ii=1:4
    avls.PRENT{ii}='d'
end
avls.tbinshft(1)=.025;
avls.tbinshft(2)=.025
avls.tbinshft(3)=.025;
avls.tbinshft(4)=.025;


avls.NFFT(1)=256;
avls.NFFT(2)=256;
avls.NFFT(3)=256;
avls.NFFT(4)=256;

avls.MKFV=[1:3]

avls.PTIND=[]
avls.SEQIND=[1:3]


avls.bastms{1}={'2007-11-19' '2007-11-20'}
avls.wntms{1}={'2007-11-21' '2007-11-27'},
avls.rectms{1}={'2007-11-28' '2007-11-30'}
avls.offset_dys=[6];
avls.LEARNANALIND=1;
% avls.bastms{2}={'2007-6-22' '2007-6-26'}
% avls.wntms{2}={'2007-6-27' '2007-7-04'}
% avls.rectms{2}={'2007-7-04' '2007-7-06'}
% 
% avls.bastms{3}={'2007-9-10' '2007-9-10'}
% avls.wntms{3}={'2007-9-11' '2007-9-16'}
% avls.rectms{3}={'2007-9-17' '2007-9-18'}
% 

avls.SEQTRGNT{1}=2;
avls.ALLSEQNT=[1 2 3 4]


