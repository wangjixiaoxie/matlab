%this script is written to generate pitch contours for bk57w35
%writes data to a .mat file in directory where data is.

%YOU NEED TO SET MUINDS BELOW
%calls jc_pitchcontourFV, in ~/matlab
%      maketimevec2
clear fv fbins bt smoothamp baspath conbins avls
avls.baspath='/cardinal3/r39g9/'

% Syntax data 
    avls.pvls{1}='probin_screen/'
    avls.cvl{1}='batchfiles'
    
    avls.WINDOWSIZE=16;
    
    avls.pvls{2}='090810_5mmap5/'
    avls.cvl{2}='batchfiles'

    avls.pvls{3}='0908_acsf/'
    avls.cvl{3}='batchfiles'

    avls.pvls{4}='0908_wnon/'
    avls.cvl{4}='batchfiles'

    avls.pvls{5}='0909_5mMap5/'
    avls.cvl{5}='batchfiles'


    avls.pvls{6}='0909_acsf/'
    avls.cvl{6}='batchfiles'

    avls.pvls{7}='0910_ap5/'
    avls.cvl{7}='batchfiles'


    avls.pvls{8}='0910_acsf/'
    avls.cvl{8}='batchfiles'

    avls.pvls{9}='0911_ap5/'
    avls.cvl{9}='batchfiles'
% 
    avls.pvls{10}='0911_acsf/'
    avls.cvl{10}='batchfiles'
    
% Which notes are hit with WN to reduce their probability?
    avls.targnts(1)='a';
    avls.targnts(2)='c';
    avls.NT{1}='a';
    avls.NT{2}='b';
    avls.NT{3}='c';
    avls.NT{4}='i'

%%%%%%%
avls.indCurve=[1 3 4 6 8 10];
avls.indCurveAC=1:length(avls.indCurve);
avls.indCurveAPV=[];
avls.indPoint=[2 5 7 9];
avls.indPointAC=[];
avls.indPointAPV=[1:4];
avls.pbase=0.7592;
