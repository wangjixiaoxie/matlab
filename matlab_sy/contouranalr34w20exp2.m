%8.19.09
clear fv fbins bt smoothamp baspath avls wnin wnrevin
%MUINDS ARE NEW DATA TO ANALYZE; PREVIOUSLY ANALYZED DATA saved to file.

avls.datdir='/oriole5/r34w20-2/datasum'
wnin(1).tmon{1}={'2009-7-26 7' }
% wnin(1).tmon{2}={'2009-7-29 7' }  
wnin(1).tmoff{1}={'2009-8-02 7'}
% wnin(1).tmoff{2}={'2009-9-04 7' '2009-9-09 7'} 

wnin(1).tmon{2}={'2009-8-09 7' }
% wnin(1).tmon{2}={'2009-7-29 7' }  
wnin(1).tmoff{2}={'2009-8-20 7'}
% wnin(1).tmoff{2}={'2009-9-04 7' '2009-9-09 7'} 


wnrevin(1).tmon{1}={'2009-8-2 7'}
wnrevin(1).tmoff{1}={'2009-8-05 7'}

wnrevin(1).tmon{2}={'2009-8-20 7'}
wnrevin(1).tmoff{2}={'2009-8-24 7'}




avls.wnin=wnin
avls.wnrevin=wnrevin

avls.pvls{1}='/oriole6/r34w20/stim_leftonly/'
avls.cvl{1}='batch'
avls.NT{1}='b'

avls.pvls{2}='/oriole5/r34w20-2/tst18_23config_90msdel_60ms_100ma/'
avls.cvl{2}='batch'
avls.NT{2}='a'


avls.pvls{3}='/oriole5/r34w20-2/tst17_100ma_60ms_13config/'
avls.cvl{3}='batch'
avls.NT{3}='a'

avls.pvls{4}='/oriole5/r34w20-2/tst7_leftonly_13_40ms_50microamp/'
avls.cvl{4}='batch'
avls.NT{4}='a'

avls.pvls{5}='/oriole5/r34w20-2/tst13_leftonly_23config//'
avls.cvl{5}='batch'
avls.NT{5}='a'


avls.pvls{6}='/oriole5/r34w20-2/tst16_leftonly_config23_100ma_40ms/'
avls.cvl{6}='batch'
avls.NT{6}='a'

avls.pvls{7}='/oriole5/r34w20-2/tst15_leftonly_13config_80ma/'
avls.cvl{7}='batch'
avls.NT{7}='a'



% avls.pvls{1}='/oriole6/r34w20/stim724_left23_right14_80ma_40ms_200Hz/'
% avls.cvl{1}='batch'
% avls.NT{1}='b'
% 
% avls.pvls{2}='/oriole6/r34w20/stim724_left23_right14_80ma_40ms_200Hz/'
% avls.cvl{2}='batch'
% avls.NT{2}='b'
% 
% 
% avls.pvls{3}='/oriole6/r34w20/bilat_14_run2_100ma_40Hz_200Hz_90msdel'
% avls.cvl{3}='batch'
% avls.NT{3}='b'
% 
% avls.pvls{4}='/oriole6/r34w20/CORRECT2_stim725_bilat14_90msdel_60ms_200Hz_60ma'
% avls.cvl{4}='batch'
% avls.NT{4}='b'
% 
% avls.pvls{5}='/oriole6/r34w20/stim725_100ma_bilat14_200Hz_40ms'
% avls.cvl{5}='batch'
% avls.NT{5}='b'
% 
% avls.pvls{6}='/oriole6/r34w20/stim729'
% avls.cvl{6}='batch'
% avls.NT{6}='b'
% 
% avls.pvls{7}='/oriole6/r34w20/stim730_90ma'
% avls.cvl{7}='batch'
% avls.NT{7}='b'
% 
% avls.pvls{8}='/oriole6/r34w20/stim731'
% avls.cvl{8}='batch'
% avls.NT{8}='b'
% 
% avls.pvls{9}='/oriole6/r34w20/stim801'
% avls.cvl{9}='batch01'
% avls.NT{9}='b'
% 
% avls.pvls{10}='/oriole6/r34w20/stim801-2'
% avls.cvl{10}='batch'
% avls.NT{10}='b'
%  
% avls.pvls{11}='/oriole6/r34w20/stimuprev_801pm'
% avls.cvl{11}='batch'
% avls.NT{11}='b'
% 
% avls.pvls{12}='/oriole6/r34w20/stim802_90ma'
% avls.cvl{12}='batch'
% avls.NT{12}='b'
%  
% avls.pvls{13}='/oriole6/r34w20/stim802_90ma'
% avls.cvl{13}='batch'
% avls.NT{13}='b'
%  
% avls.pvls{14}='/oriole6/r34w20/stim803'
% avls.cvl{14}='batch'
% avls.NT{14}='b'
% 
% avls.pvls{15}='/oriole6/r34w20/ampon803_23both'
% avls.cvl{15}='batch'
% avls.NT{15}='b'
%    
% avls.pvls{16}='/oriole6/r34w20/stim_aug17_both'
% avls.cvl{16}='batch'
% avls.NT{16}='b'
%    
% avls.pvls{17}='/oriole6/r34w20/stimon_aug18'
% avls.cvl{17}='batch'
% avls.NT{17}='b'
%  
% avls.pvls{18}='/oriole6/r34w20/ampon2_stim'
% avls.cvl{18}='batch11'
% avls.NT{18}='b'
% 
% avls.pvls{19}='/oriole6/r34w20/ampon2_stim'
% avls.cvl{19}='batch12'
% avls.NT{19}='b'
% 
% avls.pvls{20}='/oriole6/r34w20/ampon2_stim'
% avls.cvl{20}='batch14'
% avls.NT{20}='b'
% 
% avls.pvls{21}='/oriole6/r34w20/ampon2_stim'
% avls.cvl{21}='batch15'
% avls.NT{21}='b'
% 
% avls.pvls{22}='/oriole6/r34w20/stimaug19'
% avls.cvl{22}='batch'
% avls.NT{22}='b'
% 
% avls.pvls{23}='/oriole6/r34w20/stimamprevaug20'
% avls.cvl{23}='batch'
% avls.NT{23}='b'
% 
% avls.pvls{24}='/oriole6/r34w20/stimaug21'
% avls.cvl{24}='batch'
% avls.NT{24}='b'
% 
% avls.pvls{25}='/oriole6/r34w20/stim_both'
% avls.cvl{25}='batch'
% avls.NT{25}='b'
% 
% avls.pvls{26}='/oriole6/r34w20/stim_rightonly'
% avls.cvl{26}='batch'
% avls.NT{26}='b'
% 
% avls.pvls{27}='/oriole6/r34w20/stim_leftonly'
% avls.cvl{27}='batch'
% avls.NT{27}='b'

for ii=1:length(avls.pvls)
    avls.del{ii}=.090;
end
    
    
    avls.con_tbinshft=-.02;
    avls.pt_tbinshft=.015;
    avls.con_NFFT=4096;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
    avls.pt_NFFT=512;
    
    avls.conbins=[2100 3000];
    avls.fbins=[4200 6000];
%     save BINS_B NFFT fbins tbinshft
% frequency analysis just for 'b'
    avls.contanal=1:2
    avls.contnum=1:2;
    avls.contfv=[1:2]
    avls.ptfv=[1:2]
    avls.analfv=[1:2]
    avls.catchstimfv=[1:2]
    avls.mnbas=[4000]
    avls.stdbas=[100]
    avls.basruns=[];
    avls.SOUNDSTR='wav'
    avls.STIMSTR='rig'
    
    
    avls.CORSOUNDSTR='und'
    avls.CORSTIMSTR='tim'
    avls.CORSTRINDS=[]
    
    avls.STIMBNDS=[200 0]
    avls.HST_EDGES=[4200:50:6000]
    avls.STAN_RUNS=[1]
    avls.mnbas=4876
    avls.stdbas=78;
