clear fv fbins bt smoothamp baspath conbins avls wnin wnrevin

%when does wn go on and off, and asymptote??
%DELAY???  ZERO?


avls.contanal=1;
avls.datdir='/oriole6/bk57w35/datsum'

wnin(1).tmon{1}={'2009-5-20 7'}
wnin(1).tmon{2}={'2009-5-31 7'}  
wnin(1).tmoff{1}={'2009-5-24 7'}
wnin(1).tmoff{2}={'2009-6-05 7' } 
wnin(1).tmon{3}={'2009-6-21 7' } 
wnin(1).tmoff{3}={'2009-6-27 7' } 

% wnrevin(1)=[];
% wnrevin(1).tmoff{1}={}
% wnrevin(1).tmon{2}={}
% wnrevin(1).tmoff{2}={}
% wnrevin(1).tmon{3}={}
% wnrevin(1).tmoff{3}={}
avls.wnin=wnin
avls.wnrevin.tmon=[];
avls.wnrevin.tmoff=[];



avls.pvls{1}='/oriole6/bk57w35/stimoff515/'
avls.cvl{1}='batch16.catch.keep.rand'
avls.NT{1}='a'

avls.pvls{2}='/oriole6/bk57w35/520_100ma_40ms_bilat2/'
avls.cvl{2}='batch'
avls.NT{2}='a'

avls.pvls{3}='/oriole6/bk57w35/521_100ma_40ms_bilat2'
avls.cvl{3}='batch'
avls.NT{3}='a'

avls.pvls{4}='/oriole6/bk57w35/522_100ma_40ms_400Hz/'
avls.cvl{4}='batch'
avls.NT{4}='a'

avls.pvls{5}='/oriole6/bk57w35/522_100ma_40ms_400Hz_2/'
avls.cvl{5}='batch'
avls.NT{5}='a'

avls.pvls{6}='/oriole6/bk57w35/523_100ma_400Hz_40ms/'
avls.cvl{6}='batch'
avls.NT{6}='a'

avls.pvls{7}='/oriole6/bk57w35/524_100mz_100Hz_40ms/'
avls.cvl{7}='batch'
avls.NT{7}='a'

avls.pvls{8}='/oriole6/bk57w35/524_100ma_400Hz_40ms/'
avls.cvl{8}='batch'
avls.NT{8}='a'

avls.pvls{9}='/oriole6/bk57w35/526_100mz_100Hz_40ms/'
avls.cvl{9}='batch26'
avls.NT{9}='a'

avls.pvls{10}='/oriole6/bk57w35/526_100ma_400Hz_40ms/'
avls.cvl{10}='batch'
avls.NT{10}='a'


avls.pvls{11}='/oriole6/bk57w35/526_100ma_40Hz_40ms/'
avls.cvl{11}='batch'
avls.NT{11}='a'

avls.pvls{12}='/oriole6/bk57w35/528_100ma_100Hz_40ms/'
avls.cvl{12}='batch28'
avls.NT{12}='a'

avls.pvls{13}='/oriole6/bk57w35/528_100ma_400Hz_40ms/'
avls.cvl{13}='batch28'
avls.NT{13}='a'

avls.pvls{14}='/oriole6/bk57w35/529_100Hz_100ma_newconfig/'
avls.cvl{14}='batch'
avls.NT{14}='a'

avls.pvls{15}='/oriole6/bk57w35/530_newconfig_100Hz_CORRECT/'
avls.cvl{15}='batch'
avls.NT{15}='a'

avls.pvls{16}='/oriole6/bk57w35/601_400Hz_60ms_oldconfig/'
avls.cvl{16}='batch'
avls.NT{16}='a'

avls.pvls{17}='/oriole6/bk57w35/602_100Hz_60ms_origconfig/'
avls.cvl{17}='batch'
avls.NT{17}='a'

avls.pvls{18}='/oriole6/bk57w35/602_400Hz_60ms_orig_2/'
avls.cvl{18}='batch'
avls.NT{18}='a'

avls.pvls{19}='/oriole6/bk57w35/602_40Hz_60ms_origconfig/'
avls.cvl{19}='batch02'
avls.NT{19}='a'

avls.pvls{20}='/oriole6/bk57w35/603_100Hz_60ms_origconfig/'
avls.cvl{20}='batch'
avls.NT{20}='a'

avls.pvls{21}='/oriole6/bk57w35/603_100Hz_origconfig_60ms/'
avls.cvl{21}='batch'
avls.NT{21}='a'

avls.pvls{22}='/oriole6/bk57w35/603_400Hz_60ms_origconfig/'
avls.cvl{22}='batch'
avls.NT{22}='a'

avls.pvls{23}='/oriole6/bk57w35/603_40Hz_60ms_origconfig/'
avls.cvl{23}='batch'
avls.NT{23}='a'

avls.pvls{24}='/oriole6/bk57w35/617_125ma_100Hz_40ms'                                 
avls.cvl{24}='batch'
avls.NT{24}='a'

avls.pvls{25}='/oriole6/bk57w35/617_400Hz_100microamps'
avls.cvl{25}='batch'
avls.NT{25}='a'

avls.pvls{26}='/oriole6/bk57w35/617_400Hz_100microamps_40ms'               
avls.cvl{26}='batch'
avls.NT{26}='a'

avls.pvls{27}='/oriole6/bk57w35/618_innerconfig_400Hz_60ms'
avls.cvl{27}='batch'
avls.NT{27}='a'

avls.pvls{28}='/oriole6/bk57w35/618_middleconfig_400Hz_60ms'              
avls.cvl{28}='batch'
avls.NT{28}='a'

avls.pvls{29}='/oriole6/bk57w35/618_outsideconfig_400Hz_60ms'              
avls.cvl{29}='batch'
avls.NT{29}='a'

avls.pvls{30}='/oriole6/bk57w35/619_20ms_zerodel_400Hz_100ma'               
avls.cvl{30}='batch'
avls.NT{30}='a'

avls.pvls{31}='/oriole6/bk57w35/619_middlecon_400Hz_20ms_20msdel'           
avls.cvl{31}='batch'
avls.NT{31}='a'

avls.pvls{32}='/oriole6/bk57w35/619_middleconfig_100Hz_0del_60ms_12ma'    
avls.cvl{32}='batch'
avls.NT{32}='a'

avls.pvls{33}='/oriole6/bk57w35/619_middleconfig400Hz_35msdel_30ms_100ma'   
avls.cvl{33}='batch'
avls.NT{33}='a'

avls.pvls{34}='/oriole6/bk57w35/619_middleconfig_400Hz_zerodel_30ms_100ma' 
avls.cvl{34}='batch'
avls.NT{34}='a'

avls.pvls{35}='/oriole6/bk57w35/623_400Hz_midconfig_15msdel_15ms'                
avls.cvl{35}='batch'
avls.NT{35}='a'
avls.del{35}=.015;

avls.pvls{36}='/oriole6/bk57w35/623_400Hz_mid_zerodel_40ms'
avls.cvl{36}='batch'
avls.NT{36}='a'
avls.del{36}=0

avls.pvls{37}='/oriole6/bk57w35/623_outer_125ma_400ma_60ms'
avls.cvl{37}='batch'
avls.NT{37}='a'
avls.del{37}=0

avls.pvls{38}='/oriole6/bk57w35/623_outer_400Hz_60ms_zerodel'
avls.cvl{38}='batch'
avls.NT{38}='a'
avls.del{38}=0

avls.pvls{39}='/oriole6/bk57w35/623_outer_400Hz_zerodel_40ms'
avls.cvl{39}='batch'
avls.NT{39}='a'
avls.del{39}=0

avls.pvls{40}='/oriole6/bk57w35/623_variabledelays_15mspulse_400Hz_midde'
avls.cvl{40}='batch'
avls.NT{40}='a'
avls.del{40}=0

avls.pvls{41}='/oriole6/bk57w35/624_100ma_400Hz_60ms_zerodel'
avls.cvl{41}='batch'
avls.NT{41}='a'
avls.del{41}=0

avls.pvls{42}='/oriole6/bk57w35/test14'
avls.cvl{42}='batch'
avls.NT{42}='a'
avls.del{42}=0

avls.pvls{43}='/oriole6/bk57w35/test16'
avls.cvl{43}='batch'
avls.NT{43}='a'
avls.del{43}=0

avls.pvls{44}='/oriole6/bk57w35/test18'
avls.cvl{44}='batch'
avls.NT{44}='a'
avls.del{44}=0

avls.pvls{45}='/oriole6/bk57w35/623_400Hz_midconfig_25msdel_15ms'                
avls.cvl{45}='batch'
avls.NT{45}='a'
avls.del{45}=.025;

avls.pvls{46}='/oriole6/bk57w35/624var/15msdel'                
avls.cvl{46}='batch'
avls.NT{46}='a'
avls.del{46}=.015;

avls.pvls{47}='/oriole6/bk57w35/624var/20msdel'                
avls.cvl{47}='batch'
avls.NT{47}='a'
avls.del{47}=.020;

avls.pvls{48}='/oriole6/bk57w35/624var/25msdel'                
avls.cvl{48}='batch'
avls.NT{48}='a'
avls.del{48}=.025;

avls.pvls{49}='/oriole6/bk57w35/624var/30msdel'                
avls.cvl{49}='batch'
avls.NT{49}='a'
avls.del{49}=.030;

avls.pvls{50}='/oriole6/bk57w35/624var/35msdel'                
avls.cvl{50}='batch'
avls.NT{50}='a'
avls.del{50}=.035;

avls.pvls{51}='/oriole6/bk57w35/624var/40msdel'                
avls.cvl{51}='batch'
avls.NT{51}='a'
avls.del{51}=.035;



for ii=1:34
    avls.del{ii}=0
end

% clear fvpt valsa
avls.pt_NFFT=512;
avls.con_NFFT=4096;
avls.fbins=[6000 8000]
for ii=1:51
    avls.conbins{ii}=[2050 2500]
end
for ii=4:6
    avls.conbinds{ii}=[2200 2600]
end
avls.pt_tbinshft=.077;
avls.con_tbinshft=.05;
avls.contnum=3;
avls.basruns=[];

avls.contfv=[4:6]
avls.ptfv=[4:6]
avls.analfv=[1:51]
avls.catchstimfv=[51]

avls.SOUNDSTR='und'
avls.STIMSTR='im'

avls.STIMBNDS=[200 -50]
avls.HST_EDGES=[6300:50:8000]
avls.STAN_RUNS=[4 5 6 8 10 13 16 22 24:39 41 ]

avls.mnbas=6889
avls.stdbas=104.3




