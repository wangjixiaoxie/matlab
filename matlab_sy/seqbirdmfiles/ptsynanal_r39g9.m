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
avls.baspath='/oriole7/dir2/r39g9/'


%START OFF BY LABELING KEY DATA

%pre data=11/20
pvl{1}='screen'
cvl{1}='batch14.catch'

pvl{2}='wnon_716'
cvl{2}='batch18.catch'

% pvl{3}='firstshift/seqshift_02-13-09'
% cvl{3}='batch.rnd'
% 
% pvl{4}='firstshift/seqshift_02-14-09'
% cvl{4}='batch.rnd'
% 
% pvl{5}='firstshift/seqshift_02-15-09'
% cvl{5}='batch.rnd'
% 
% pvl{6}='firstshift/seqshift_02-16-09'
% cvl{6}='batch.rnd'
% 
% pvl{7}='firstshift/seqshift_02-17-09'
% cvl{7}='batchnotes'
% 
% pvl{8}='firstshift/seqshift_02-18-09'
% cvl{8}='batchnotes'
% 
% pvl{9}='firstshift/seqshift_02-19-09'
% cvl{9}='batch.rnd'
% 
% pvl{10}='firstshift/seqshift_02-20-09'
% cvl{10}='batchnotes'
% 
% pvl{11}='firstshift/seqshift_02-21-09'
% cvl{11}='batchnotes'
% 
% pvl{12}='firstshift/seqshift_02-22-09'
% cvl{12}='batchnotes'
% 
% pvl{13}='firstshift/seqshift_02-23-09'
% cvl{13}='batchnotesc'
% 
% pvl{14}='firstshift/seqshift_02-24-09'
% cvl{14}='batchnotescbin'
% 
% pvl{15}='firstshift/recov_02-25-09'
% cvl{15}='batch.rnd'
% 
% pvl{16}='firstshift/recov_02-26-09'
% cvl{16}='batch.rnd'
% 
% pvl{17}='firstshift/recov_02-27-09'
% cvl{17}='batch.rnd'
% 
% pvl{18}='firstshift/recov_02-28-09'
% cvl{18}='batch.rnd'
% 
% pvl{19}='firstshift/recov_03-01-09'
% cvl{19}='batch.rnd'
% 
% pvl{20}='secondshift/pre_050309'
% cvl{20}='batchnotes'
% 
% 
% pvl{21}='secondshift/shift_050509'
% cvl{21}='batchnotes'
% 
% pvl{22}='secondshift/shift_050609'
% cvl{22}='batchnotes'
% 
% pvl{23}='secondshift/shift_051109'
% cvl{23}='batchnotes'
% 
% 
% 
% 
% 

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
% avls.fbins{4}=[2000 3000];

avls.NT{1}='a';
avls.NT{2}='b';
avls.NT{3}='i';
% avls.NT{3}='z'
% avls.NT{4}='i'
% 
for ii=1:3
    avls.PRENT{ii}=''
end
avls.tbinshft(1)=.025;
avls.tbinshft(2)=.025
avls.tbinshft(3)=.025;
% avls.tbinshft(4)=.025;


avls.NFFT(1)=256;
avls.NFFT(2)=256;
avls.NFFT(3)=256;
% avls.NFFT(4)=256;

avls.MKFV=[1:2]

avls.PTIND=[1:2]
avls.SEQIND=[1:2]


avls.bastms{1}={'2009-2-10' '2009-2-12'}
avls.wntms{1}={'2009-2-13' '2009-2-24'},
avls.rectms{1}={'2009-2-25' '2009-3-1'}
avls.offset_dys=[];

avls.bastms{2}={'2009-5-01' '2009-5-03'}
avls.wntms{2}={'2009-5-04' '2009-5-14'},
avls.rectms{2}={'2009-5-15' '2009-5-20'}
avls.offset_dys=[];


% avls.bastms{2}={'2007-6-22' '2007-6-26'}
% avls.wntms{2}={'2007-6-27' '2007-7-04'}
% avls.rectms{2}={'2007-7-04' '2007-7-06'}
% 
% avls.bastms{3}={'2007-9-10' '2007-9-10'}
% avls.wntms{3}={'2007-9-11' '2007-9-16'}
% avls.rectms{3}={'2007-9-17' '2007-9-18'}
% 

avls.SEQTRGNT=1;
avls.ALLSEQNT=[1 2 ]


