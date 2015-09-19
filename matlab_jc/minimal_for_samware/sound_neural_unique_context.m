%
%
%
function sig_out=sound_neural_unique_context(file_number,syllable_to_quant,premotorwind)


if strfind(pwd,'pu77bk41')
    birdname='pu77bk41';
elseif strfind(pwd,'pu24w39')
    birdname='pu24w39';
elseif strfind(pwd,'pu26y2')
    birdname='pu26y2';
elseif strfind(pwd,'bl82bl81')
    birdname='bl82bl81';
elseif strfind(pwd,'B53O71')
    birdname='B53O71';
elseif strfind(pwd,'r93bl81')
    birdname='r93bl81';
elseif strfind(pwd,'o85pu54')
    birdname='o85pu54';
elseif strfind(pwd,'G45G46')
    birdname='G45G46';
elseif strfind(pwd,'Pu55Pu22')
    birdname='Pu55Pu22';
elseif strfind(pwd,'W90O73')
    birdname='W90O73';
elseif strfind(pwd,'W32Pi51')
    birdname='W32Pi51';
elseif strfind(pwd,'W15W94')
    birdname='W15W94';
elseif strfind(pwd,'O14O15')
    birdname='O14O15';
elseif strfind(pwd,'Pk35G27')
    birdname='Pk35G27';
elseif strfind(pwd,'W96Pi45')
    birdname='W96Pi45';
elseif strfind(pwd,'g38o18') %mel lman
    birdname='g38o18';
elseif strfind(pwd,'p85g54')%mel lman
    birdname='p85g54';
elseif strfind(pwd,'b39b14')%mimi lman
    birdname='b39b14';
elseif strfind(pwd,'b9r63')%mimi lman
    birdname='b9r63';
elseif strfind(pwd,'p49w84')%mimi lman
    birdname='p49w84';
elseif strfind(pwd,'blueblue')%mimi lman
    birdname='blueblue';
elseif strfind(pwd,'jman3')%mimi lman
    birdname='jman3';
elseif strfind(pwd,'g91pu54')%mel lman
    birdname='g91pu54';
elseif strfind(pwd,'masa')%mimi lman
    birdname='masa';
elseif strfind(pwd,'pur98w3')%mimi lman
    birdname='pur98w3';
elseif strfind(pwd,'pu77bk41')%mimi lman
    birdname='pu77bk41';
elseif strfind(pwd,'o85bk54')%mimi lman
    birdname='o85bk54';
    
end

clear global spikes_on_x_axis
global spikes_on_x_axis
spikes_on_x_axis=1;

%premotor_window=[30 70];
%premotor_window=[20 60];
premotor_window=premotorwind;  % NOTE:  This window gives significant correlations when looking at HVC spike data and BOS amplitude.  
%premotor_window=[60 100];disp('USING ACAUSAL WINDOW')
[seq_array,skip_first,skip_last,t_assay]=generate_seq_array(birdname,syllable_to_quant);

analysis_type_code_B53O71=2;
analysis_type_code_r93bl81=2;
analysis_type_code_o85pu54=2;
analysis_type_code_o85pu54=2;
analysis_type_code_G45G46=2;
analysis_type_code_Pu55Pu22=2;
analysis_type_code_W90O73=2;
analysis_type_code_W32Pi51=2;
analysis_type_code_W15W94=2;
analysis_type_code_O14O15=2;
analysis_type_code_Pk35G27=2;
analysis_type_code_W96Pi45=2;
analysis_type_code_g38o18=2;
analysis_type_code_p85g54=2;
analysis_type_code_b39b14=2;
analysis_type_code_p85g54=2;
analysis_type_code_g91pu54=2;
analysis_type_code_b39b14=2;
analysis_type_code_p49w84=2;
analysis_type_code_b9r63=2;
analysis_type_code_blueblue=2;
analysis_type_code_jman3=2;
analysis_type_code_masa=2;
analysis_type_code_pu77bk41=2;
analysis_type_code_o85bk54=2;

pur98w3_array{1}='combined_data_pur98w3_SU26_PCA_CH1_TH_recommended.mat';
pur98w3_array{end+1}='combined_data_pur98w3_SU29_PCA_CH0_TH_recommended.mat';

masa_array{1}='combined_data_masa_051700_PCA_CH1_TH_recommended.mat';
masa_array{end+1}='combined_data_masa_072000_PCA_CH1_TH_recommended.mat';
masa_array{end+1}='combined_data_masa_073100_PCA_CH0_TH_recommended.mat';

jman3_array{1}='combined_data_jman3_061802_PCA_CH1_TH_recommended.mat';
jman3_array{end+1}='combined_data_jman3_070102_PCA_CH1_TH_recommended.mat';

blueblue_array{1}='combined_data_blueblue_site1_PCA_CH0_TH_recommended.mat';

b9r63_array{1}='combined_data_b9r63_3_PCA_CH1_TH_recommended.mat';

p49w84_array{1}='combined_data_p49w84_SU5_PCA_CH0_TH_recommended.mat';
p49w84_array{end+1}='combined_data_p49w84_SU6_PCA_CH0_TH_recommended.mat';
p49w84_array{end+1}='combined_data_p49w84_SU7_PCA_CH0_TH_recommended.mat';

b39b14_array{1}='combined_data_b39b14_SU11_PCA_CH0_TH_recommended.mat';
b39b14_array{end+1}='combined_data_b39b14_SU12_PCA_CH0_TH_recommended.mat';
b39b14_array{end+1}='combined_data_b39b14_SU16_PCA_CH0_TH_recommended.mat';


g38o18_array{1}='combined_data_g38o18_mu1_e1_PCA_CH0_TH_recommended.mat';
g38o18_array{end+1}='combined_data_g38o18_mu2_e2_PCA_CH1_TH_recommended.mat';
g38o18_array{end+1}='combined_data_g38o18_mu3_e1_PCA_CH0_TH_recommended.mat';
g38o18_array{end+1}='combined_data_g38o18_mu4_e1_PCA_CH0_TH_recommended.mat';
g38o18_array{end+1}='combined_data_g38o18_su1_e1_PCA_CH0_TH_recommended.mat';
g38o18_array{end+1}='combined_data_g38o18_su2_e1_PCA_CH0_TH_recommended.mat';

g91pu54_array{1}='combined_data_g91pu54unit1_PCA_CH0_TH_recommended.mat';
% g91pu54_array{end+1}='combined_data_g91pu54_mu2_e1_PCA_CH0_TH_recommended.mat';
% g91pu54_array{end+1}='combined_data_g91pu54_mu3_e1_PCA_CH0_TH_recommended.mat';
% g91pu54_array{end+1}='combined_data_g91pu54_su1_e2_PCA_CH1_TH_recommended.mat';

p85g54_array{1}='combined_data_p85g54_mu1_e2_PCA_CH1_TH_recommended.mat';
p85g54_array{end+1}='combined_data_p85g54_mu2_e2_PCA_CH0_TH_recommended.mat';

if 1% PCA-based techniaue
    %     analysis_type_code=1;  Older threshold technique
    %     analysis_type_code=2;  Newer PCA-based technique
    %
    analysis_type_code_pu77bk41=2;
    pu77bk41_array{1}='combined_data_pu77bk41unit1_PCA_CH0_TH_recommended.mat';
%     pu77bk41_array{end+1}='combined_data_pu77bk41_3_17_2006_1409-1544_PCA_CH1_TH_recommended.mat';
%     pu77bk41_array{end+1}='combined_data_pu77bk41_3_17_2006_1600_1614_PCA_CH1_TH_recommended.mat';
%     pu77bk41_array{end+1}='combined_data_pu77bk41_3_17_2006_1600_1614_PCA_CH2_TH_recommended.mat';
%     pu77bk41_array{end+1}='combined_data_pu77bk41_3_17_2006_1600_1614_PCA_CH3_TH_recommended.mat';
%     pu77bk41_array{end+1}='combined_data_pu77bk41_3_17_2006_1805-2025_PCA_CH2_TH_reco.mat';
%     pu77bk41_array{end+1}='combined_data_pu77bk41_3_17_2006_1805-2025_PCA_CH3_TH_rec.mat';
%     pu77bk41_array{end+1}='combined_data_pu77bk41_3_19_2005_0808-0839_PCA_CH1_TH_recommended.mat';
%     pu77bk41_array{end+1}='combined_data_pu77bk41_3_19_2006_1752-1831_PCA_CH3_TH_recommended.mat';
%     pu77bk41_array{end+1}='combined_data_pu77bk41_3_19_2006_1752-1831_PCA_CH4_TH_recommended.mat';
%     pu77bk41_array{end+1}='combined_data_pu77bk41_3_20_2006_1026_1125_PCA_CH2_TH_recommended.mat';
%     pu77bk41_array{end+1}='combined_data_pu77bk41_3_20_2006_1622-1632_PCA_CH2_TH_recommended.mat';
%     pu77bk41_array{end+1}='combined_data_pu77bk41_3_20_2006_1634-1645_PCA_CH2_TH_recommended.mat';
%     pu77bk41_array{end+1}='combined_data_pu77bk41_3_20_3006_1344-1413_PCA_CH2_TH_recommended.mat';
%     pu77bk41_array{end+1}='combined_data_pu77bk41_3_20_3006_1516-1521_PCA_CH1_TH_recommended.mat';
%     pu77bk41_array{end+1}='combined_data_pu77bk41_3_21_2006_1321-1340_PCA_CH1_TH_recommended.mat';
%     pu77bk41_array{end+1}='combined_data_pu77bk41_3_21_2006_1321-1340_PCA_CH2_subset_TH_recommended.mat';
%     pu77bk41_array{end+1}='combined_data_pu77bk41_3_21_2006_1702-1745_PCA_CH1_TH_recommended.mat';
%     pu77bk41_array{end+1}='combined_data_pu77bk41_3_22_2006_0945-0958_PCA_CH1_TH_recommended.mat';
%     pu77bk41_array{end+1}='combined_data_pu77bk41_3_22_2006_1014-1018_PCA_CH3_TH_recommended.mat';
%     pu77bk41_array{end+1}='combined_data_pu77bk41_3_22_2006_1020-1029_PCA_CH2_TH_recommended.mat';
%     pu77bk41_array{end+1}='combined_data_pu77bk41_3_22_2006_1020-1029_PCA_CH4_TH_recommended.mat';
%     pu77bk41_array{end+1}='combined_data_pu77bk41_3_24_2006_1428-1443_PCA_CH4_TH_recommended.mat';
%     pu77bk41_array{end+1}='combined_data_pu77bk41_3_24_2006_1448-1454_PCA_CH1_TH_recommended.mat';
%     pu77bk41_array{end+1}='combined_data_pu77bk41_3_24_2006_1448-1454_PCA_CH3_TH_recommended.mat';
else% old (threshold) way
    analysis_type_code_pu77bk41=1;
    pu77bk41_array{1}='combined_data_pu77bk41_3_17_2006_1805-2025';
    pu77bk41_array{end+1}='combined_data_pu77bk41_3_20_3006_1344-1413_ch2';
    pu77bk41_array{end+1}='combined_data_pu77bk41_3_20_3006_1516-1521_ch1';
    pu77bk41_array{end+1}= 'combined_data_pu77bk41_3_17_2006_1409-1544';
    pu77bk41_array{end+1}= 'combined_data_pu77bk41_3_20_2006_1026_1125_ch2';
    pu77bk41_array{end+1}='combined_data_pu77bk41_3_17_2006_2025-2042_multi_ch2';
    pu77bk41_array{end+1}='combined_data_pu77bk41_3_22_2006_0945-0958_bigMulti_ch1';
    %% below are multiunits
    pu77bk41_array{end+1}= 'combined_data_pu77bk41_3_17_2006_1600_1614_multi_ch1';
    pu77bk41_array{end+1}=  'combined_data_pu77bk41_3_17_2006_1600_1614_multi_ch2';
    pu77bk41_array{end+1}= 'combined_data_pu77bk41_3_17_2006_1600_1614_multi_ch3';
    pu77bk41_array{end+1}= 'combined_data_pu77bk41_3_17_2006_1805-2025_multi_ch3';
    pu77bk41_array{end+1}=  'combined_data_pu77bk41_3_19_2005_0808-0839_multi_ch1';
    pu77bk41_array{end+1}= 'combined_data_pu77bk41_3_19_2006_1752-1831_multi_ch3';
    pu77bk41_array{end+1}=  'combined_data_pu77bk41_3_19_2006_1752-1831_multi_ch4';
    pu77bk41_array{end+1}=  'combined_data_pu77bk41_3_20_2006_1622-1632_multi_ch2';
    pu77bk41_array{end+1}=  'combined_data_pu77bk41_3_20_2006_1634-1645_bigmulti_ch2';
    pu77bk41_array{end+1}= 'combined_data_pu77bk41_3_21_2006_1321-1340_multi_ch1';
    pu77bk41_array{end+1}=  'combined_data_pu77bk41_3_21_2006_1321-1340_multi_ch2';
    pu77bk41_array{end+1}='combined_data_pu77bk41_3_21_2006_1702-1745_multi_ch2';
    pu77bk41_array{end+1}='combined_data_pu77bk41_3_22_2006_1014-1018_multi_ch3';
    pu77bk41_array{end+1}='combined_data_pu77bk41_3_22_2006_1020-1029_multi_ch2';
    pu77bk41_array{end+1}='combined_data_pu77bk41_3_22_2006_1020-1029_multi_ch4';
    pu77bk41_array{end+1}='combined_data_pu77bk41_3_24_2006_1428-1443_multi_ch4';
    pu77bk41_array{end+1}='combined_data_pu77bk41_3_24_2006_1448-1454_multi_ch1';
    pu77bk41_array{end+1}='combined_data_pu77bk41_3_24_2006_1448-1454_multi_ch3';
end

analysis_type_code_bl82bl81=2;
bl82bl81_array{1}='combined_data_bl82bl81unit1_PCA_CH1_TH_recommended.mat';
% bl82bl81_array{1}='combined_data_bl82bl81unit2_PCA_CH1_TH_recommended.mat';
% bl82bl81_array{end+1}='combined_data_bl82bl81_2_15_2006_1012-1038_PCA_CH1_TH_recommended.mat';
% bl82bl81_array{end+1}='combined_data_bl82bl81_2_5_2006_1352-1437_PCA_CH1_TH_recommended.mat';
% bl82bl81_array{end+1}='combined_data_bl82bl81_2_5_2006_1437-1510_PCA_CH1_TH_recommended.mat';
% bl82bl81_array{end+1}='combined_data_bl82bl81_2_5_2006_1437-1510_PCA_CH2_TH_recommended.mat';
% bl82bl81_array{end+1}='combined_data_bl82bl81_2_5_2006_1437-1510_PCA_CH4_TH_recommended.mat';
% bl82bl81_array{end+1}='combined_data_bl82bl81_2_5_2006_1525-1635_PCA_CH1_TH_recommended.mat';
% bl82bl81_array{end+1}='combined_data_bl82bl81_2_5_2006_1525-1635_PCA_CH4_TH_recommended.mat ';
% bl82bl81_array{end+1}='combined_data_bl82bl81_2_5_2006_1635-0847_PCA_CH4_TH_recommended.mat';
% bl82bl81_array{end+1}='combined_data_bl82bl81_2_5_2006_1635-0847_PCA_CH3_TH_recommended.mat';
% bl82bl81_array{end+1}='combined_data_bl82bl81_2_6_2006_0852-0910_PCA_CH3_TH_recommended.mat';
% bl82bl81_array{end+1}='combined_data_bl82bl81_2_6_2006_0923-1018_PCA_CH1_TH_recommended.mat';



Pu55Pu22_array{1}='combined_data_Pu55Pu22_MU_01_19_2007_1935_1954_PCA_CH0_TH_recommended.mat';
Pu55Pu22_array{end+1}='combined_data_Pu55Pu22_MU_01_19_2007_1955_2126_PCA_CH0_TH_recommended.mat';
Pu55Pu22_array{end+1}='combined_data_Pu55Pu22_MU_01_20_2007_1325_1359_PCA_CH0_TH_recommended.mat';
Pu55Pu22_array{end+1}='combined_data_Pu55Pu22_MU_01_20_2007_1838_2048_PCA_CH0_TH_recommended.mat';
Pu55Pu22_array{end+1}='combined_data_Pu55Pu22_MU_01_20_2007_1442_1523_PCA_CH0_TH_recommended.mat';
Pu55Pu22_array{end+1}='combined_data_Pu55Pu22_MU_01_21_2007_1352_2019_PCA_CH0_TH_recommended.mat';
Pu55Pu22_array{end+1}='combined_data_Pu55Pu22_MU_01_22_2007_1235_1355_PCA_CH0_TH_recommended.mat';
Pu55Pu22_array{end+1}='combined_data_Pu55Pu22_MU_01_22_2007_1523_1702_PCA_CH0_TH_recommended.mat';
Pu55Pu22_array{end+1}='combined_data_Pu55Pu22_MU_01_22_2007_1828_2012_PCA_CH0_TH_recommended.mat';
Pu55Pu22_array{end+1}='combined_data_Pu55Pu22_SU_01_22_2007_1235_1632_PCA_CH1_TH_recommended.mat';
Pu55Pu22_array{end+1}='combined_data_Pu55Pu22_MU_01_23_2007_1139_1300_PCA_CH0_TH_recommended.mat';
Pu55Pu22_array{end+1}='combined_data_Pu55Pu22_MU_01_23_2007_1336_1430_PCA_CH0_TH_recommended.mat';
Pu55Pu22_array{end+1}='combined_data_Pu55Pu22_MU_01_23_2007_1336_1430_PCA_CH2_TH_recommended.mat';
Pu55Pu22_array{end+1}='combined_data_Pu55Pu22_MU_01_23_2007_1551_1859_PCA_CH2_TH_recommended.mat';
%Pu55Pu22_array{end+1}='combined_data_Pu55Pu22_MU_01_19_2007_1250_1440_PCA_
%CH0_TH_recommended.mat';  % IN-like


B53O71_array{1}='combined_data_PCA__CH0_TH_recommended.mat';
B53O71_array{end+1}='combined_data_B53O71_MU_05_05_2006_1000_1057_PCA_CH0_TH_recommended.mat';


o85pu54_array{1}='combined_data_o85pu54unit1_PCA_CH0_TH_recommended.mat';
% o85pu54_array{end+1}='combined_data_o85pu54_MU_5_30_2006_1034_1113_PCA_CH0_TH_recommended.mat';
% o85pu54_array{end+1}='combined_data_o85pu54_SU_5_28_2006_1020_1314_PCA_CH0_TH_recommended.mat';
% o85pu54_array{end+1}='combined_data_o85pu54_MU_5_29_2006_1530_1540_PCA_CH0_TH_recommended.mat';
%o85pu54_array{end+1}='combined_data_o85pu54_SU_5_30_2006_1211_1258_PCA_CH0_T
%H_recommended.mat'; % INTERNEURON-LIKE
% below is newly added
% o85pu54_array{end+1}='combined_data_o85pu54_MU_5_30_2006_0734_0941_PCA_CH0_TH_recommended.mat';



W90O73_array{1}='combined_data_W90O73_MU_10_11_2005_1851_2147_PCA_CH0_TH_recommended.mat';
W90O73_array{end+1}='combined_data_W90O73_MU_10_1213_2005_0949_0957_PCA_CH0_TH_recommended.mat';

G45G46_array{1}='combined_data_G45G46_MU_03_06_2004_0952_1049_PCA_CH1_TH_recommended.mat';


W32Pi51_array{1}='combined_data_W32Pi51_MU_10_09_2006_1656_2100_PCA_CH2_TH_recommended.mat';
W32Pi51_array{end+1}='combined_data_W32Pi51_MU_10_09_2006_1656_2100_PCA_CH3_TH_recommended.mat';
W32Pi51_array{end+1}='combined_data_W32Pi51_MU_10_09_2006_0703_0725_PCA_CH2_TH_recommended.mat';
W32Pi51_array{end+1}='combined_data_W32Pi51_MU_10_09_2006_0728_0906_PCA_CH2_TH_recommended.mat';
W32Pi51_array{end+1}='combined_data_W32Pi51_MU_10_09_2006_0728_0906_PCA_CH3_TH_recommended.mat';
W32Pi51_array{end+1}='combined_data_W32Pi51_MU_10_09_2006_0936_1013_PCA_CH2_TH_recommended.mat';
W32Pi51_array{end+1}='combined_data_W32Pi51_MU_10_10_2006_1109_1137_PCA_CH2_TH_recommended.mat';
W32Pi51_array{end+1}='combined_data_W32Pi51_MU_10_10_2006_1109_1137_PCA_CH3_TH_recommended.mat';
W32Pi51_array{end+1}='combined_data_W32Pi51_MU_10_10_2006_1145_1304_PCA_CH2_TH_recommended.mat';
W32Pi51_array{end+1}='combined_data_W32Pi51_MU_10_10_2006_1320_1357_PCA_CH2_TH_recommended.mat';
W32Pi51_array{end+1}='combined_data_W32Pi51_MU_10_10_2006_1401_1452_PCA_CH2_TH_recommended.mat';
W32Pi51_array{end+1}='combined_data_W32Pi51_MU_10_10_2006_1401_1452_PCA_CH3_TH_recommended.mat';
W32Pi51_array{end+1}='combined_data_W32Pi51_MU_10_10_2006_1952_2022_PCA_CH2_TH_recommended.mat';
W32Pi51_array{end+1}='combined_data_W32Pi51_MU_10_10_2006_1952_2022_PCA_CH3_TH_recommended.mat';


r93bl81_array{1}='combined_data_r93bl81unit1_PCA_CH0_TH_recommended.mat';
% r93bl81_array{end+1}='combined_data_r93bl81_MU_08_12_2006_1010_1044_PCA_CH3_TH_recommended.mat';
% r93bl81_array{end+1}='combined_data_r93bl81_MU_08_13_2006_0714_0718_PCA_CH0_TH_recommended.mat';
% r93bl81_array{end+1}='combined_data_r93bl81_MU_08_13_2006_0702_0719_PCA_CH2_TH_recommended.mat';
% r93bl81_array{end+1}='combined_data_r93bl81_SU_08_11_2006_0900_0930_PCA_CH0_TH_recommended.mat';
% r93bl81_array{end+1}='combined_data_r93bl81_SU_08_11_2006_0702_0708_PCA_CH0_TH_recommended.mat';
% r93bl81_array{end+1}='combined_data_r93bl81_SU_08_11_2006_1805_1837_PCA_CH0_TH_recommended.mat';
% r93bl81_array{end+1}='combined_data_r93bl81_SU_08_11_2006_1805_1837_PCA_CH3_TH_recommended.mat';

W15W94_array{1}='combined_data_W15W94_MU_08_08_2003_1401_9999_PCA_CH0_TH_recommended.mat';
W15W94_array{end+1}='combined_data_W15W94_MU_12_01_2003_1122_1614_PCA_CH0_TH_recommended.mat';


O14O15_array{1}='combined_data_O14O15_MU_02_17_2006_1451_1547_PCA_CH0_TH_recommended.mat';
O14O15_array{end+1}='combined_data_O14O15_MU_02_18_2006_1112_1204_PCA_CH0_TH_recommended.mat';
O14O15_array{end+1}='combined_data_O14O15_MU_02_19_2006_1354_1541_PCA_CH0_TH_recommended.mat';
O14O15_array{end+1}='combined_data_O14O15_MU_02_27_2006_1350_1858_PCA_CH1_TH_recommended.mat';
O14O15_array{end+1}='combined_data_O14O15_MU_02_28_2006_1253_1611_PCA_CH1_TH_recommended.mat';
O14O15_array{end+1}='combined_data_O14O15_MU_03_01_2006_1346_2125_PCA_CH1_TH_recommended.mat';
O14O15_array{end+1}='combined_data_O14O15_MU_03_02_2006_1150_1305_PCA_CH1_TH_recommended.mat';


Pk35G27_array{1}='combined_data_Pk35G27_MU_05_29_2007_1750_1757_PCA_CH0_TH_recommended.mat';
Pk35G27_array{end+1}='combined_data_Pk35G27_MU_05_29_2007_1750_1757_PCA_CH1_TH_recommended.mat';
Pk35G27_array{end+1}='combined_data_Pk35G27_MU_05_29_2007_1750_1757_PCA_CH2_TH_recommended.mat';
Pk35G27_array{end+1}='combined_data_Pk35G27_MU_05_30_2007_0714_0748_PCA_CH0_TH_recommended.mat';
Pk35G27_array{end+1}='combined_data_Pk35G27_MU_05_30_2007_0714_0748_PCA_CH1_TH_recommended.mat';
Pk35G27_array{end+1}='combined_data_Pk35G27_MU_05_30_2007_0714_0748_PCA_CH2_TH_recommended.mat';
Pk35G27_array{end+1}='combined_data_Pk35G27_MU_05_30_2007_0714_0748_PCA_CH3_TH_recommended.mat';
Pk35G27_array{end+1}='combined_data_Pk35G27_MU_05_30_2007_0823_1029_PCA_CH1_TH_recommended.mat';
Pk35G27_array{end+1}='combined_data_Pk35G27_MU_05_30_2007_0823_1029_PCA_CH3_TH_recommended.mat';
Pk35G27_array{end+1}='combined_data_Pk35G27_MU_05_30_2007_1220_1608_PCA_CH0_TH_recommended.mat';
Pk35G27_array{end+1}='combined_data_Pk35G27_MU_05_30_2007_1220_1608_PCA_CH1_TH_recommended.mat';
Pk35G27_array{end+1}='combined_data_Pk35G27_MU_05_30_2007_1220_1608_PCA_CH3_TH_recommended.mat';
Pk35G27_array{end+1}='combined_data_Pk35G27_MU_05_30_2007_1519_1523_PCA_CH3_TH_recommended.mat';
Pk35G27_array{end+1}='combined_data_Pk35G27_MU_05_31_2007_0716_0925_PCA_CH1_TH_recommended.mat';
Pk35G27_array{end+1}='combined_data_Pk35G27_MU_05_31_2007_0716_0925_PCA_CH2_TH_recommended.mat';
Pk35G27_array{end+1}='combined_data_Pk35G27_MU_05_31_2007_0959_1040_PCA_CH0_TH_recommended.mat';
Pk35G27_array{end+1}='combined_data_Pk35G27_MU_05_31_2007_1144_1426_PCA_CH0_TH_recommended.mat';
Pk35G27_array{end+1}='combined_data_Pk35G27_MU_05_31_2007_1144_1426_PCA_CH1_TH_recommended.mat';
Pk35G27_array{end+1}='combined_data_Pk35G27_MU_06_01_2007_1130_1215_PCA_CH0_TH_recommended.mat';
Pk35G27_array{end+1}='combined_data_Pk35G27_MU_06_01_2007_1130_1215_PCA_CH1_TH_recommended.mat';
Pk35G27_array{end+1}='combined_data_Pk35G27_MU_06_01_2007_1130_1215_PCA_CH2_TH_recommended.mat';
Pk35G27_array{end+1}='combined_data_Pk35G27_MU_06_01_2007_1237_1344_PCA_CH1_TH_recommended.mat';

W96Pi45_array{1}='combined_data_W96Pi45_MU_06_03_2007_2021_0858_PCA_CH0_TH_recommended.mat';
W96Pi45_array{end+1}='combined_data_W96Pi45_MU_06_04_2007_1254_0804_PCA_CH2_TH_recommended.mat';
W96Pi45_array{end+1}='combined_data_W96Pi45_MU_06_05_2007_0925_1135_PCA_CH2_TH_recommended.mat';
W96Pi45_array{end+1}='combined_data_W96Pi45_MU_06_05_2007_1147_1304_PCA_CH3_TH_recommended.mat';
W96Pi45_array{end+1}='combined_data_W96Pi45_MU_06_05_2007_1320_1539_PCA_CH1_TH_recommended.mat';
W96Pi45_array{end+1}='combined_data_W96Pi45_MU_06_05_2007_1320_1539_PCA_CH3_TH_recommended.mat';
%W96Pi45_array{end+1}='combined_data_W96Pi45_MU_06_05_2007_1659_1754_PCA_CH
%1_TH_recommended.mat';  % same as a SU
W96Pi45_array{end+1}='combined_data_W96Pi45_MU_06_05_2007_1844_2027_PCA_CH3_TH_recommended.mat';
W96Pi45_array{end+1}='combined_data_W96Pi45_MU_06_13_2007_1401_2057_PCA_CH1_TH_recommended.mat';
W96Pi45_array{end+1}='combined_data_W96Pi45_MU_06_13_2007_1401_2057_PCA_CH2_TH_recommended.mat';
W96Pi45_array{end+1}='combined_data_W96Pi45_MU_06_13_2007_1401_2057_PCA_CH3_TH_recommended.mat';
W96Pi45_array{end+1}='combined_data_W96Pi45_MU_06_14_2007_0702_0728_PCA_CH0_TH_recommended.mat';
W96Pi45_array{end+1}='combined_data_W96Pi45_MU_06_14_2007_1707_1826_PCA_CH1_TH_recommended.mat';
W96Pi45_array{end+1}='combined_data_W96Pi45_MU_06_14_2007_1707_1826_PCA_CH2_TH_recommended.mat';
W96Pi45_array{end+1}='combined_data_W96Pi45_SU_06_05_2007_1755_1808_PCA_CH1_TH_recommended.mat';
W96Pi45_array{end+1}='combined_data_W96Pi45_SU_06_13_2007_1336_1358_PCA_CH2_TH_recommended.mat';

% o85pu54_array{1}='combined_data_o85pu54_MELDATA_5_28_2006_SU';
% o85pu54_array{end+1}='combined_data_o85pu54_MELDATA_5_30_2006_SU_CH3r';
% o85pu54_array{end+1}='combined_data_o85pu54_MELDATA_5_30_2006_SU_CH4r';
% o85pu54_array{end+1}='combined_data_o85pu54_MELDATA_5_28_2006_MU_INCOMPLETE';
% o85pu54_array{end+1}='combined_data_o85pu54_MELDATA_5_29_2006_SU';
% o85pu54_array{end+1}='combined_data_o85pu54_MELDATA_5_29_2006_MU_f1_to_f29';

if 1
    if 1
        analysis_type_code_pu26y2=2;
        pu26y2_array{1}='combined_data_pu26y2_11_13_2006_1341_1654_PCA_CH1_TH_recommended.mat';
        pu26y2_array{end+1}='combined_data_pu26y2_11_13_2006_1743_1854_PCA_CH1_TH_recommended.mat';
        pu26y2_array{end+1}='combined_data_pu26y2_11_17_2006_1023_1155_PCA_CH1_TH_recommended.mat';
        pu26y2_array{end+1}='combined_data_pu26y2_11_20_2006_1019_1212_PCA_CH1_TH_recommended.mat';
        pu26y2_array{end+1}='combined_data_pu26y2_11_22_2006_1431_1826_PCA_CH1_TH_recommended.mat';
        pu26y2_array{end+1}='combined_data_pu26y2_11_30_2006_1522_9999_PCA_CH1_TH_recommended.mat';
        pu26y2_array{end+1}='combined_data_pu26y2_12_01_2006_1054_1107_PCA_CH1_TH_recommended.mat';
        pu26y2_array{end+1}='combined_data_pu26y2_12_01_2006_1242_1305_PCA_CH1_TH_recommended.mat';
        pu26y2_array{end+1}='combined_data_pu26y2_12_01_2006_1306_1330_PCA_CH1_TH_recommended.mat';
        pu26y2_array{end+1}='combined_data_pu26y2_12_01_2006_1401_1415_PCA_CH1_TH_recommended.mat';
        pu26y2_array{end+1}='combined_data_pu26y2_12_01_2006_1415_1559_PCA_CH1_TH_recommended.mat';
    elseif 0    % second side -RA?
        analysis_type_code_pu26y2=2;
        pu26y2_array{1}='combined_data_pu26y2_12_27_2006_1352_1429_PCA_CH3_TH_recommended.mat';
        pu26y2_array{end+1}='combined_data_pu26y2_01_04_2007_1342_1439_PCA_CH4_TH_recommended.mat';
        pu26y2_array{end+1}='combined_data_pu26y2_01_10_2007_1511_1541_PCA_CH1_TH_recommended.mat';
        pu26y2_array{end+1}='combined_data_pu26y2_01_10_2007_1559_1623_PCA_CH4_TH_recommended.mat';
        pu26y2_array{end+1}='combined_data_pu26y2_01_12_2007_1413_1443_PCA_CH4_TH_recommended.mat';
        pu26y2_array{end+1}='combined_data_pu26y2_01_12_2007_1456_1528_PCA_CH4_TH_recommended.mat';
        pu26y2_array{end+1}='combined_data_pu26y2_01_16_2007_1058_1214_PCA_CH1_TH_recommended.mat';
        % pu26y2_array{end+1}='
        % pu26y2_array{end+1}='
        % pu26y2_array{end+1}='
    end
else
    analysis_type_code_pu26y2=1;
    % all recordins are from ch1
    % after line - %[#songs #syl j]
    pu26y2_array{1}='combined_data_pu26y2_11_13_2006_1341_1654_undirected';     %[13  76]
    pu26y2_array{end+1}='combined_data_pu26y2_11_13_2006_1743_1854_undirected'; %[14  71]
    pu26y2_array{end+1}='combined_data_pu26y2_11_17_2006_1023_1155_directed';   %[6   33]
    pu26y2_array{end+1}='combined_data_pu26y2_11_17_2006_1023_1155_undirected'; %[44 239]
    %pu26y2_array{end+1}='combined_data_pu26y2_11_20_2006_1019_1212_directed';   %[1    4]
    pu26y2_array{end+1}='combined_data_pu26y2_11_20_2006_1019_1212_undirected'; %[17 105]

    %pu26y2_array{end+1}='combined_data_pu26y2_11_22_2006_1021_1026_undirected'; %[4   19]
    %pu26y2_array{end+1}='combined_data_pu26y2_11_22_2006_1026_1100_undirected'; %[2   13]
    %pu26y2_array{end+1}='combined_data_pu26y2_11_22_2006_1431_1826_directed';   %[2    5]
    pu26y2_array{end+1}='combined_data_pu26y2_11_22_2006_1431_1826_undirected'; %[33 218]
    %pu26y2_array{end+1}='combined_data_pu26y2_11_28_2006_1437_1505_directed';   %[2   14]
    %pu26y2_array{end+1}='combined_data_pu26y2_11_28_2006_1437_1505_undirected'; %[2   10]
    %pu26y2_array{end+1}='combined_data_pu26y2_11_28_2006_1534_1604_directed';   %[1    4]
    %pu26y2_array{end+1}='combined_data_pu26y2_11_28_2006_1534_1604_undirected'; %[2   10]
    %pu26y2_array{end+1}='combined_data_pu26y2_11_30_2006_1522_9999_directed';   %[1   12]
    pu26y2_array{end+1}='combined_data_pu26y2_11_30_2006_1522_9999_undirected'; %[11  75]
    pu26y2_array{end+1}='combined_data_pu26y2_12_01_2006_1415_1559_undirected'; %[21 166]
    %pu26y2_array{end+1}='combined_data_pu26y2_12_01_2006_1057_1306_undirected';
    pu26y2_array{end+1}='combined_data_pu26y2_12_01_2006_1242_1305_partial_undirected';
    pu26y2_array{end+1}='combined_data_pu26y2_12_01_2006_1054_1107_partial_undirected';

    pu26y2_array{end+1}='combined_data_pu26y2_12_01_2006_1306_1330_undirected';
    pu26y2_array{end+1}='combined_data_pu26y2_12_01_2006_1401_1415_undirected';

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%55
    % second implant
    pu26y2_array{end+1}='combined_data_pu26y2_12_27_2006_1352_1429_undirected'; % SU from ch3 - not as sparse as usual
    pu26y2_array{end+1}='combined_data_pu26y2_01_04_2007_1342_1439_undirected'; %might not be motor - one burst, but sustained tonic
    pu26y2_array{end+1}='combined_data_pu26y2_01_10_2007_1511_1541_undirected';  % this and below seem a bit atypical for RA
    pu26y2_array{end+1}='combined_data_pu26y2_01_10_2007_1559_1623_ch4_undirected';

    pu26y2_array{end+1}='combined_data_pu26y2_01_12_2007_1413_1443_undirected';
    pu26y2_array{end+1}='combined_data_pu26y2_01_12_2007_1456_1528_undirected';

    pu26y2_array{end+1}='combined_data_pu26y2_01_16_2007_1058_1214_undirected';  % again, not typical RA burstiness.  nice single unit though
end

eval(sprintf('analysis_type_code=analysis_type_code_%s',birdname))
if ~isstr(file_number)
    if strcmp(birdname,'pu77bk41')
        loadfilename=pu77bk41_array{file_number};
    elseif strcmp(birdname,'o85pu54') % (mel's bird)
        loadfilename=o85pu54_array{file_number};
    elseif strcmp(birdname,'pu26y2')
        loadfilename=pu26y2_array{file_number};
    elseif strcmp(birdname,'bl82bl81')% NEW
        loadfilename=bl82bl81_array{file_number};
    end
    disp(['Loading ' loadfilename])%return
    sig_out=analyze(seq_array,skip_first,skip_last,t_assay,loadfilename,spikes_on_x_axis,premotor_window);
    %        set(gcf,'position',[-3 158 1280 804]) % fullscreen
elseif strcmp(file_number,'all')        % do for all files and save
    eval(sprintf('use_array=%s_array;',birdname))
    ct=1;
    ofst=700;
    for x=1:length(use_array)
        figure(ct+ofst);clf;
        disp(['Loading ' use_array{x}])
        sig_out(ct,:)=analyze(seq_array,skip_first,skip_last,t_assay,use_array{x},spikes_on_x_axis,premotor_window,syllable_to_quant);
        ct=ct+1;
    end
    if 1    % save output
        eval(sprintf('sig_out_%s = sig_out',syllable_to_quant))
        if strcmp(syllable_to_quant,'B')
            % Matlab wont overwrite a filename differing only in Caps.  So
            % for 'B' vs 'b', save 'B' as 'Bcap'
            eval(sprintf('save sig_out_%s_syl_Bcap sig_out_%s use_array premotor_window analysis_type_code',birdname,syllable_to_quant))
        else
            %            eval(sprintf('save sig_out_%s_syl_%s sig_out_%s use_array premotor_window analysis_type_code',birdname,syllable_to_quant,syllable_to_quant))
            eval(sprintf('save sig_out_%s_syl_%s_20_60_msec sig_out_%s use_array premotor_window analysis_type_code',birdname,syllable_to_quant,syllable_to_quant))
        end
    else
        disp('  ');disp(' NOT SAVING OUTPUT ');disp('  ')
    end
else
    disp(['Loading ' file_number])%return
    sig_out=analyze(seq_array,skip_first,skip_last,t_assay,file_number,spikes_on_x_axis,premotor_window,syllable_to_quant);
end

file_number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function sig_out=analyze(seq_array,skip_first,skip_last,t_assay,loadfilename,spikes_on_x_axis,premotor_window,syllable_to_quant)

n_factors=4;% just pitch, amp

eval(sprintf('load %s',loadfilename))
%sprintf('load %s_low_pitch_c',loadfilename)
%eval(sprintf('load %s_low_pitch_c',loadfilename))
ct=1;
for y=1:length(seq_array)
    clear true
    for x=1:(length(labels)-length(seq_array{y})+1)
        compare=labels(x:(x+length(seq_array{y})-1));
        true(x)=strcmp(seq_array{y},compare);
    end
    n_true=length(find(true));
    %    disp('  ')
    if n_true>0
        disp(' USING ANY SEQ > 1 ')
    end
    %   disp('  ')

    if n_true>=1
        %        if n_true==15
        use_seq_array{ct}=seq_array{y};
        ct=ct+1;
    end
    disp(['Found ' num2str(n_true)  ' instances of ' seq_array{y}])
end


symb_vec='osdv^v<>ph';
symb_vec=[symb_vec symb_vec];

% initialize empty variables
all_spikes=[];all_pitches=[];all_amplitudes=[];all_entropies=[];all_harmonic_ratios=[];all_serial_orders=[];all_durations=[];
all_amplitude_slopes=[];all_pitch_slopes=[];all_entropy_slopes=[];

all_spikes_minusmean=[];all_pitches_minusmean=[];all_amplitudes_minusmean=[];all_entropies_minusmean=[];
all_harmonic_ratios_minusmean=[];all_serial_orders_minusmean=[];all_durations_minusmean=[];all_amplitude_slopes_minusmean=[];
all_pitch_slopes_minusmean=[];all_entropy_slopes_minusmean=[];

mean_spikes=[];mean_pitch=[];mean_amplitudes=[];mean_entropies=[];mean_harmonic_ratios=[];
mean_serial_orders=[];mean_durations=[];mean_amplitude_slopes=[];mean_pitch_slopes=[];mean_entropy_slopes=[];

for x=1:length(use_seq_array)

    %    [spikes_tmp,pitch_tmp,amplitude_tmp,entropy_tmp,harmonics_tmp,serial_orders_tmp,mean_outs]=foo2(use_seq_array{x},skip_first,skip_last,loadfilename,symb_vec(x),t_assay,premotor_window);
    [spikes_data,behavioral_data,mean_outs]=foo2(use_seq_array{x},skip_first,skip_last,loadfilename,symb_vec(x),t_assay,premotor_window);

    % contents of spikes_data and behavioral_data for seq_array of length = y
    %     x_data{y}=spikes
    %     y_data{y,1}=pitch
    %     y_data{y,2}=amplitude
    %     y_data{y,3}=entropy
    %     y_data{y,4}=serial_order
    %     y_data{y,5}=harmonics
    %     y_data{y,6}=duration
    %     y_data{y,7}=amplitude slope
    %     y_data{y,8}=pitch slope
    %     y_data{y,9}=entropy slope

    tmp=size(behavioral_data,1);
    for z=1:tmp
        %        harmonics_tmp{z}=log10(harmonics_tmp{z});
        harmonics_tmp{z}=log10(behavioral_data{z,5});

        if x*z==1;disp('TAKING LOG OF HARMONICS RATIO');end
    end
    for y=1:length(use_seq_array{x})-skip_first-skip_last

        spikes_minusmean=spikes_data{y}-mean(spikes_data{y});
        pitch_minusmean=behavioral_data{y,1}-mean(behavioral_data{y,1});
        amplitude_minusmean=behavioral_data{y,2}-mean(behavioral_data{y,2});
        entropy_minusmean=behavioral_data{y,3}-mean(behavioral_data{y,3});
        serial_orders_minusmean=behavioral_data{y,4}-mean(behavioral_data{y,4});
        harmonic_ratios_minusmean=behavioral_data{y,5}-mean(behavioral_data{y,5});
        durations_minusmean=behavioral_data{y,6}-mean(behavioral_data{y,6});
        amplitude_slopes_minusmean=behavioral_data{y,7}-mean(behavioral_data{y,7});
%         pitch_slopes_minusmean=behavioral_data{y,8}-mean(behavioral_data{y,8});
%         entropy_slopes_minusmean=behavioral_data{y,9}-mean(behavioral_data{y,9});

        all_spikes=[all_spikes spikes_data{y}'];
        all_pitches=[all_pitches behavioral_data{y,1}'];
        all_amplitudes=[all_amplitudes behavioral_data{y,2}'];
        all_entropies=[all_entropies behavioral_data{y,3}'];
        all_serial_orders=[all_serial_orders behavioral_data{y,4}'];
        all_harmonic_ratios=[all_harmonic_ratios behavioral_data{y,5}'];
        all_durations=[all_durations behavioral_data{y,6}'];
        all_amplitude_slopes=[all_amplitude_slopes behavioral_data{y,7}'];
%         all_pitch_slopes=[all_pitch_slopes behavioral_data{y,8}'];
%         all_entropy_slopes=[all_entropy_slopes behavioral_data{y,9}'];

        all_spikes_minusmean=[all_spikes_minusmean spikes_minusmean'];
        all_pitches_minusmean=[all_pitches_minusmean pitch_minusmean'];
        all_amplitudes_minusmean=[all_amplitudes_minusmean amplitude_minusmean'];
        all_entropies_minusmean=[all_entropies_minusmean entropy_minusmean'];
        all_harmonic_ratios_minusmean=[all_harmonic_ratios_minusmean harmonic_ratios_minusmean'];
        all_serial_orders_minusmean=[all_serial_orders_minusmean serial_orders_minusmean'];
        all_durations_minusmean=[all_durations_minusmean durations_minusmean'];
        all_amplitude_slopes_minusmean=[all_amplitude_slopes_minusmean amplitude_slopes_minusmean'];
%         all_pitch_slopes_minusmean=[all_pitch_slopes_minusmean pitch_slopes_minusmean'];
%         all_entropy_slopes_minusmean=[all_entropy_slopes_minusmean entropy_slopes_minusmean'];

        mean_spikes=[mean_spikes mean_outs.spikes];
        mean_pitch=[mean_pitch mean_outs.pitch];
        mean_amplitudes=[mean_amplitudes mean_outs.amp];
        mean_entropies=[mean_entropies mean_outs.entropy];
        mean_harmonic_ratios=[mean_harmonic_ratios mean_outs.harmonics];
        mean_serial_orders=[mean_serial_orders mean_outs.serial_order];
        mean_durations=[mean_durations mean_outs.duration];
        mean_amplitude_slopes=[mean_amplitude_slopes mean_outs.amp_slope];
%         mean_pitch_slopes=[mean_pitch_slopes mean_outs.pitch_slope];
%         mean_entropy_slopes=[mean_entropy_slopes mean_outs.entropy_slope];
    end
end

length_all_pitches=length(all_pitches)

subplot(3,4,1);hold on;
%subplot(3,3,1);hold on;
%subplot(3,2,1);hold on;
if 0
    if spikes_on_x_axis
        plot_and_regress(mean_spikes',mean_pitch',[1 1 1]*.5,1,[],[]);
    else
        plot_and_regress(mean_pitch',mean_spikes',[1 1 1]*.5,1,[],[]);
    end
end

subplot(3,4,2);hold on;
%subplot(3,3,2);hold on;
%subplot(3,2,2);hold on;
if 0
    if spikes_on_x_axis
        plot_and_regress(mean_spikes',mean_amplitudes',[1 1 1]*.5,1,[],[]);
    else
        plot_and_regress(mean_amplitudes',mean_spikes',[1 1 1]*.5,1,[],[]);
    end
end
clear xlt
xlt{1}=RemoveUnderScore(loadfilename);
xlabel(xlt)
%t_str{1}=['All correlations computed in window onset + [' num2str(pre_onset) '  '   num2str(post_onset) '] msec' ];
t_str{1}=['All correlations computed relative to t\_assay = ' num2str(t_assay) ];
if length(use_seq_array)==1
    t_str{2}=['Looking for sequence ' use_seq_array{1} ', skipping first ' num2str(skip_first) ' and last ' num2str(skip_last) ];
end
title(t_str)

%%%%%%%%%%%%%%%%%%%%%   DURATION %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%   NOT plotted %%%%%%%%%%%%%%%%%
if 0
    tmpfig=gcf;
    figure(100);clf;hold on

    x_data=all_spikes';y_data=all_durations';

    [tmp,tmp2,regress_struct_duration]=plot_and_regress(x_data,y_data,[1 0 0],1,'# spikes','duration');
    figure(tmpfig)
    %regress_struct_duration,return
else
    if spikes_on_x_axis
        x_data=all_spikes';y_data=all_durations';
        % remove outliers based on y data
        [x_data,y_data]=remove_outliers(x_data,y_data);
    else
        y_data=all_durations';x_data=all_entropies';
    end

    [FIT,tmp,tmp,tmp,STATS] =regress(y_data,[x_data ones(size(x_data))]);
    regress_struct_duration.FIT=FIT;regress_struct_duration.STATS=STATS;
    %    regress_struct_duration,return
end

%%%%%%%%%%%%%%%%%%%%%   AMPLITUDE SLOPE %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%   NOT plotted %%%%%%%%%%%%%%%%%
%     tmpfig=gcf;figure(100);clf;hold on;x_data=all_spikes';y_data=all_amplitude_slopes';
%     [tmp,tmp2,regress_struct_duration]=plot_and_regress(x_data,y_data,[1 0 0],1,'# spikes','amp slope'); figure(tmpfig)
if spikes_on_x_axis
    x_data=all_spikes';y_data=all_amplitude_slopes';
    % remove outliers based on y data
    [x_data,y_data]=remove_outliers(x_data,y_data);
else
    y_data=all_amplitude_slopes';x_data=all_entropies';
end

[FIT,tmp,tmp,tmp,STATS] =regress(y_data,[x_data ones(size(x_data))]);
regress_struct_all_amplitude_slope.FIT=FIT;regress_struct_all_amplitude_slope.STATS=STATS;

% %%%%%%%%%%%%%%%%%%%%%   PITCH SLOPE %%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%   NOT plotted %%%%%%%%%%%%%%%%%
% if spikes_on_x_axis
%     x_data=all_spikes';y_data=all_pitch_slopes';
%     % remove outliers based on y data
%     [x_data,y_data]=remove_outliers(x_data,y_data);
% else
%     y_data=all_pitch_slopes';x_data=all_entropies';
% end
% 
% [FIT,tmp,tmp,tmp,STATS] =regress(y_data,[x_data ones(size(x_data))]);
% regress_struct_all_pitch_slope.FIT=FIT;regress_struct_all_pitch_slope.STATS=STATS;
% 
% %%%%%%%%%%%%%%%%%%%%%   ENTROPY SLOPE %%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%   NOT plotted %%%%%%%%%%%%%%%%%
% if spikes_on_x_axis
%     x_data=all_spikes';y_data=all_entropy_slopes';
%     % remove outliers based on y data
%     [x_data,y_data]=remove_outliers(x_data,y_data);
% else
%     y_data=all_entropy_slopes';x_data=all_entropies';
% end
% 
% [FIT,tmp,tmp,tmp,STATS] =regress(y_data,[x_data ones(size(x_data))]);
% regress_struct_all_entropy_slope.FIT=FIT;regress_struct_all_entropy_slope.STATS=STATS;

%     tmpfig=gcf;figure(100);clf;hold on;
%     [tmp,tmp2,tmp3]=plot_and_regress(x_data,y_data,[1 0 0],1,'# spikes','entropy slope'); figure(tmpfig)


if n_factors>2
    %%%%%%%%%%%%%%%%%%%%%   ENTROPY %%%%%%%%%%%%%%%%%
    subplot(3,4,7);cla;hold on;
    %subplot(3,3,6);cla;hold on;
    if spikes_on_x_axis
        [tmp,tmp2,regress_struct_entropy]=plot_and_regress(all_spikes',all_entropies',[1 0 0],1,'# spikes','spectral entropy');
    else
        [tmp,tmp2,regress_struct_entropy]=plot_and_regress(all_pitches',all_entropies',[1 0 0],1,'spectral entropy','# spikes');
    end

    subplot(3,4,11);hold on;
    %subplot(3,3,9);hold on;
    if spikes_on_x_axis
        plot_and_regress(all_spikes_minusmean',all_entropies_minusmean',[0 0 1],1,'# spikes, minus within-context mean','entropy, minus within-context mean');
    else
        plot_and_regress(all_entropies_minusmean',all_spikes_minusmean',[0 0 1],1,'entropy, minus within-context mean','# spikes, minus within-context mean');
    end
end

if n_factors>3
    %%%%%%%%%%%%%%%%%%%%%   HARMONIC POWER %%%%%%%%%%%%%%%%%
    subplot(3,4,8);cla;hold on;
    if spikes_on_x_axis
        %        [tmp,tmp2,regress_struct_hi_lo_harmonics]=plot_and_regress(all_spikes',all_harmonic_ratios',[1 0 0],1,'# spikes','Power ratio hi/lo harmonic');
        [tmp,tmp2,regress_struct_serial_order]=plot_and_regress(all_spikes',all_serial_orders',[1 0 0],1,'# spikes','Serial order');

        if ~sum(abs(rem(all_serial_orders,1)))  % if all serial orders are integers
            for y=unique(all_serial_orders)
                id=find(all_serial_orders==y);
                if length(id)>1
                    plot(mean(all_spikes(id)),y,'b+','linew',2)
                elseif 0 % to plot context text
                    if length(seq_array)==1 & length(seq_array{1})==1
                        id2=find(labels==seq_array{1});
                        text(all_spikes(id)+1,y,labels(id2(y)-2:id2(y)+2))
                    end
                end
            end
        end

    else
        %        [tmp,tmp2,regress_struct_hi_lo_harmonics]=plot_and_regress(all_harmonic_ratios',all_spikes',[1 0 0],1,'Power ratio hi/lo harmonic','# spikes');
        [tmp,tmp2,regress_struct_serial_order]=plot_and_regress(all_serial_orders',all_spikes',[1 0 0],1,'Serial order','# spikes');
    end

    subplot(3,4,12);hold on;
    if spikes_on_x_axis
        %        plot_and_regress(all_spikes_minusmean',all_harmonic_ratios_minusmean',[0 0 1],1,'# spikes, minus within-context mean','Power ratio hi/lo harmonic, minus within-context mean');
        plot_and_regress(all_spikes_minusmean',all_serial_orders_minusmean',[0 0 1],1,'# spikes, minus within-context mean','Serial order');
    else
        %        plot_and_regress(all_harmonic_ratios_minusmean',all_spikes_minusmean',[0 0 1],1,'Power ratio hi/lo harmonic, minus within-context mean','# spikes, minus within-context mean');
        plot_and_regress(all_serial_orders_minusmean',all_spikes_minusmean',[0 0 1],1,'Serial order, minus within-context mean','# spikes, minus within-context mean');
    end
end
%%%%%%%%%%%%%%%%%%%%%   PITCH %%%%%%%%%%%%%%%%%

subplot(3,4,5);cla;hold on;
%subplot(3,3,4);cla;hold on;
%subplot(3,2,3);cla;hold on;
if spikes_on_x_axis
    [tmp,tmp2,regress_struct]=plot_and_regress(all_spikes',all_pitches',[1 0 0],1,'# spikes','pitch (hz)');
else
    [tmp,tmp2,regress_struct]=plot_and_regress(all_pitches',all_spikes',[1 0 0],1,'pitch (hz)','# spikes');
end
% save temp_sound_neural_pitch_spikes all_pitches all_spikes
% [rho,p]=corr(all_pitches',all_spikes');
% xlabel(num2str([rho p]))
% disp('HACK')

% tmp=gcf;
% figure(100);clf;hold on;
% subplot(1,2,1);hold on
% hist(all_spikes,0:9);xlabel('# spikes in premotor window');ylabel('N')
% set(gca,'xtick',[0 5 10],'xlim',[-1 11])
% axis square
% subplot(1,2,2);hold on
% hist(all_pitches/2,2900:50:3700);xlabel('Freq (Hz)');ylabel('N')
% set(gca,'xtick',[3000:100:3600],'xlim',[2900 3700])
% axis square
% figure(tmp)

subplot(3,4,9);hold on;
%subplot(3,3,7);hold on;
%subplot(3,2,5);hold on;
if spikes_on_x_axis
    plot_and_regress(all_spikes_minusmean',all_pitches_minusmean',[0 0 1],1,'# spikes, minus within-context mean','pitch, minus within-context mean');
else
    plot_and_regress(all_pitches_minusmean',all_spikes_minusmean',[0 0 1],1,'pitch, minus within-context mean','# spikes, minus within-context mean');
end


%%%%%%%%%%%%%%%%%%%%%   AMPLITUDE %%%%%%%%%%%%%%%%%

subplot(3,4,6);hold on;
%subplot(3,3,5);hold on;
%subplot(3,2,4);hold on;
if spikes_on_x_axis
    [tmp,tmp2,regress_struct_amp]=plot_and_regress(all_spikes',all_amplitudes',[1 0 0],1,'# spikes','amplitude (a.u.)');
else
    [tmp,tmp2,regress_struct_amp]=plot_and_regress(all_amplitudes',all_spikes',[1 0 0],1,'amplitude','# spikes');
end

subplot(3,4,10);hold on;
%subplot(3,3,8);hold on;
%subplot(3,2,6);hold on;
if spikes_on_x_axis
    plot_and_regress(all_spikes_minusmean',all_amplitudes_minusmean',[0 0 1],1,'# spikes, minus within-context mean','amplitude, minus within-context mean');
else
    plot_and_regress(all_amplitudes_minusmean',all_spikes_minusmean',[0 0 1],1,'amplitude, minus within-context mean','# spikes, minus within-context mean');
end


%  output significance level, positive if positive slope, negative
% if negative slope
if regress_struct_amp.STATS(3)==0
    regress_struct_amp.STATS(3)=1e-15;
    disp('Calculated p<0.  Reseting p=1e-15')
end
if regress_struct.STATS(3)==0
    regress_struct.STATS(3)=1e-15;
    disp('Calculated p<0.  Reseting p=1e-15')
end
sig_out(1,1)=sign(regress_struct.FIT(1))*regress_struct.STATS(3);  % sign of slope times p level - pitch
sig_out(1,2)=sign(regress_struct_amp.FIT(1))*regress_struct_amp.STATS(3);  % sign of slope times p level - amplitude
if exist('SNR')
    disp('Columns 3 and 7 are SNR')
    sig_out(1,3)=mean(SNR);% mean sig to noise
    sig_out(1,7)=length(SNR);% length of songs
else
    disp('Columns 3 and 7 are pct_error')
    sig_out(1,3)=mean(pct_error);% mean sig to noise
    sig_out(1,7)=length(pct_error);% length of songs
end
sig_out(1,4)=mean(all_spikes)/((abs(diff(premotor_window)))*.001);% mean spike rate
sig_out(1,5)=regress_struct.STATS(1);% Rsquared  - pitch
sig_out(1,6)=regress_struct_amp.STATS(1);% Rsquared - amp
sig_out(1,8)=length(all_spikes);% length of syllables
sig_out(1,9)=regress_struct.FIT(1);% slope of fit to pitch
sig_out(1,10)=regress_struct_amp.FIT(1);% slope of fit to amp

if n_factors>2   %| strcmp(file_number,'all')
    sig_out(1,11)=sign(regress_struct_entropy.FIT(1))*regress_struct_entropy.STATS(3);  % sign of slope times p level - entropy
    sig_out(1,12)=regress_struct_entropy.STATS(1);% Rsquared  - entropy
    sig_out(1,13)=regress_struct_entropy.FIT(1);% slope of fit to entropy
end

if n_factors>3 %| strcmp(file_number,'all')
    sig_out(1,14)=sign(regress_struct_serial_order.FIT(1))*regress_struct_serial_order.STATS(3);  % sign of slope times p level - serial order
    sig_out(1,15)=regress_struct_serial_order.STATS(1);% Rsquared  - serial order
    sig_out(1,16)=regress_struct_serial_order.FIT(1);% slope of fit to serial order

    %     sig_out(1,14)=sign(regress_struct_hi_lo_harmonics.FIT(1))*regress_struct_hi_lo_harmonics.STATS(3);  % sign of slope times p level - harmonics ratio
    %     sig_out(1,15)=regress_struct_hi_lo_harmonics.STATS(1);% Rsquared  - harmonics ratio
    %     sig_out(1,16)=regress_struct_hi_lo_harmonics.FIT(1);% slope of fit to harmonics ratio
end

sig_out(1,17)=sign(regress_struct_duration.FIT(1))*regress_struct_duration.STATS(3);  % sign of slope times p level - duration
sig_out(1,18)=regress_struct_duration.STATS(1);% Rsquared  - duration
sig_out(1,19)=regress_struct_duration.FIT(1);% slope of fit to duration

sig_out(1,20)=sign(regress_struct_all_amplitude_slope.FIT(1))*regress_struct_all_amplitude_slope.STATS(3);  % sign of slope times p level - amplitude slope
sig_out(1,21)=regress_struct_all_amplitude_slope.STATS(1);% Rsquared  - amplitude slope
sig_out(1,22)=regress_struct_all_amplitude_slope.FIT(1);% slope of fit to amplitude slope

% sig_out(1,23)=sign(regress_struct_all_pitch_slope.FIT(1))*regress_struct_all_pitch_slope.STATS(3);  % sign of slope times p level - pitch slope
% sig_out(1,24)=regress_struct_all_pitch_slope.STATS(1);% Rsquared  - pitch slope
% sig_out(1,25)=regress_struct_all_pitch_slope.FIT(1);% slope of fit to pitch slope
% 
% sig_out(1,26)=sign(regress_struct_all_entropy_slope.FIT(1))*regress_struct_all_entropy_slope.STATS(3);  % sign of slope times p level - entropy slope
% sig_out(1,27)=regress_struct_all_entropy_slope.STATS(1);% Rsquared  - entropy slope
% sig_out(1,28)=regress_struct_all_entropy_slope.FIT(1);% slope of fit to entropy slope

id1=strfind(loadfilename,'_PCA');
id2=strfind(loadfilename,'_CH');chnum=loadfilename(id2+3);
save_name=['all_scalars_' loadfilename(15:id1-1) '_skipping' '_CH' chnum '_syl_' syllable_to_quant];
%save_name=['all_scalars_' loadfilename(15:id1-1) '_CH' chnum '_syl_' syllable_to_quant];
mean_error=mean(pct_error);% mean sig to noise
n_songs=length(pct_error);% length of songs
mean_rate=mean(all_spikes)/((abs(diff(premotor_window)))*.001);% mean spike rate

eval(sprintf('save %s all_spikes all_pitches all_amplitudes all_entropies all_pitch_slopes all_amplitude_slopes all_entropy_slopes mean_error n_songs premotor_window mean_rate' ,save_name))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%function [spikes_output,pitch_output,amplitude_output,entropy_output,harmonics_output,serial_order_output,mean_outs]=foo2(search_for,skip_first,skip_last,loadfilename,symbol,t_assay,premotor_window)
function [x_data,y_data,mean_outs]=foo2(search_for,skip_first,skip_last,loadfilename,symbol,t_assay,premotor_window)
eval(sprintf('load %s',loadfilename))


global spikes_on_x_axis
premotor_estimate=premotor_window;

spiketimes=round(spiketimes*1000);   % in milliseconds
spiketimes=unique(spiketimes);  % if more than one spike in a bin, call it one spike

sp=zeros(1,length(spiketimes));

spiketimes=nonzeros(spiketimes);

sp(spiketimes)=1000;

onsets=unique(round(onsets));
offsets=unique(round(offsets));

for x=1:(length(labels)-length(search_for)+1)
    compare=labels(x:(x+length(search_for)-1));
    true(x)=strcmp(search_for,compare);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get spike number and syl duration
id_true=find(true);                                             % id of starts of sequence
%id_true,search_for,return
disp(['Found ' num2str(length(id_true))  ' sequences'])
MAT_duration=zeros(length(id_true),length(search_for));
MAT_int_duration=MAT_duration;
MAT_nspikes=MAT_duration;
MAT_nspikes_int=MAT_duration;
for x=1:length(id_true)
    for y=1:length(search_for)

        % for calculating syllable durations
        syllable_onsets=onsets(id_true(x)+y-1);
        syllable_offsets=offsets(id_true(x)+y-1);
        if x<length(id_true)
            syllable_onsets_nextsyl=onsets(id_true(x)+y);
        else
            syllable_onsets_nextsyl=offsets(id_true(x)+y-1);
        end

        % this is where window for spike counting is set
        %disp(['At window setting - [pre_onset post_onset] = ' num2str([pre_onset post_onset])])
        % OLD  - this was wrong - pre_ and post_onset should be relative to
        % time of quantification, not to onset of syllable
        %         on=onsets(id_true(x)+y-1)-pre_onset;
        %         off=onsets(id_true(x)+y-1)+post_onset;

        if t_assay<1
            pre_onset=-1*(t_assay*1000-premotor_estimate(2));
            post_onset=t_assay*1000-premotor_estimate(1);
        else   % if t_assay>1, it is the percent time through the syllable.
            t_assay_with_pct=.01*t_assay*[offsets(id_true(x)+y-1)-onsets(id_true(x)+y-1)];
            pre_onset=-1*(t_assay_with_pct-premotor_estimate(2));
            post_onset=t_assay_with_pct-premotor_estimate(1);
            % %assuming that avg duration is 50ms.
            %             pre_onset=-1*(.01*t_assay*50-premotor_estimate(2));
            %             post_onset=.01*t_assay*50-premotor_estimate(1);
        end

        % this is where window for spike counting is set
        on=onsets(id_true(x)+y-1)-pre_onset;
        off=onsets(id_true(x)+y-1)+post_onset;

        if numel(onsets)>=(id_true(x)+y-1)
            MAT_nspikes(x,y)=length(find(spiketimes>on & spiketimes<off));   % number of spikes between on and off;

            % old measures when on and off corresponded to onset and offset
            % of syl
            MAT_duration(x,y)=(syllable_offsets-syllable_onsets);
            MAT_int_duration(x,y)=(syllable_onsets_nextsyl-syllable_offsets);

        end
    end
    MAT_amplitude(x,:)=peak_amplitudes(id_true(x):(id_true(x)+length(search_for)-1));
    MAT_amplitude_summed(x,:)=summed_amplitudes(id_true(x):(id_true(x)+length(search_for)-1));
    if exist('amplitude_at_pitchquant')
        MAT_amplitude_at_pitchquant(x,:)=amplitude_at_pitchquant(id_true(x):(id_true(x)+length(search_for)-1));
    end
    MAT_amplitude_16msec_around_pitchquant(x,:)=amplitude_16msec_around_pitchquant(id_true(x):(id_true(x)+length(search_for)-1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get pitch
for y=(1+skip_first):(length(search_for)-skip_last)
    for x=1:length(id_true)
        pitch(x,y)=peak_pinterp_labelvec(id_true(x)+y-1);
        pitch_raw(x,y)=peak_labelvec(id_true(x)+y-1);
        if x<length(id_true)
            pitch_next_syl(x,y)=peak_pinterp_labelvec(id_true(x)+y);
        else
            pitch_next_syl(x,y)=0;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get entropy
for y=(1+skip_first):(length(search_for)-skip_last)
    if exist('spectral_entropy')
        for x=1:length(id_true)
            entropy(x,y)=spectral_entropy(id_true(x)+y-1);
        end
    else
        entropy(x,y)=zeros(size( pitch_next_syl(x,y)));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get hi/lo harmonics ratio
for y=(1+skip_first):(length(search_for)-skip_last)
    if exist('hi_to_lo_harmonic_ratio')
        for x=1:length(id_true)
            harmonics(x,y)=hi_to_lo_harmonic_ratio(id_true(x)+y-1);
        end
    else
        harmonics(x,y)=zeros(size( pitch_next_syl(x,y)));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get serial_order
for y=(1+skip_first):(length(search_for)-skip_last)
    if exist('serial_order')
        for x=1:length(id_true)
            serial_order_tmp(x,y)=serial_order(id_true(x)+y-1);
        end
    else
        serial_order_tmp(x,y)=zeros(size( pitch_next_syl(x,y)));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get duration
for y=(1+skip_first):(length(search_for)-skip_last)
    for x=1:length(id_true)
        duration(x,y)=offsets(id_true(x)+y-1)-onsets(id_true(x)+y-1);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get amp_slope
for y=(1+skip_first):(length(search_for)-skip_last)
    if exist('amplitude_slope_at_pitchquant')
        for x=1:length(id_true)
            amp_slope(x,y)=amplitude_slope_at_pitchquant(id_true(x)+y-1);
        end
    else
        amp_slope(x,y)=zeros(size( pitch_next_syl(x,y)));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get pitch_slope
% for y=(1+skip_first):(length(search_for)-skip_last)
%     if exist('pitch_slope')
%         for x=1:length(id_true)
%             pitch_slope_tmp(x,y)=pitch_slope(id_true(x)+y-1);
%         end
%     else
%         pitch_slope_tmp(x,y)=zeros(size( pitch_next_syl(x,y)));
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get entropy_slope
% for y=(1+skip_first):(length(search_for)-skip_last)
%     if exist('entropy_slope')
%         for x=1:length(id_true)
%             entropy_slope_tmp(x,y)=entropy_slope(id_true(x)+y-1);
%         end
%     else
%         entropy_slope_tmp(x,y)=zeros(size( pitch_next_syl(x,y)));
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sutract mean of that syllable
for y=(1+skip_first):(length(search_for)-skip_last)

    MAT_amplitude_minusmean(:,y)=MAT_amplitude(:,y)-mean(MAT_amplitude(:,y));
    MAT_amplitude_summed_minusmean(:,y)=MAT_amplitude_summed(:,y)-mean(MAT_amplitude_summed(:,y));
    pitch_minusmean(:,y)=pitch(:,y)-mean(pitch(:,y));
    pitch_next_syl_minusmean(:,y)=pitch_next_syl(:,y)-mean(pitch_next_syl(:,y));

    MAT_nspikes_minusmean(:,y)=MAT_nspikes(:,y)-mean(MAT_nspikes(:,y));
    %    MAT_nspikes_int_prev_minusmean(:,y)=MAT_nspikes_int_prev(:,y)-mean(MAT_nspikes_int_prev(:,y));
    pitch_minusmean(:,y)=pitch(:,y)-mean(pitch(:,y));

    MAT_duration_minusmean(:,y)=MAT_duration(:,y)-mean(MAT_duration(:,y));
    MAT_int_duration_minusmean(:,y)=MAT_int_duration(:,y)-mean(MAT_int_duration(:,y));

    MAT_serial_pos(:,y)=y*ones(length(MAT_nspikes),1);

end

if find(MAT_duration==0);    error('zero in MAT_duration');end

MAT_mean_amplitude=MAT_amplitude./MAT_duration;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% Plotting
syl_vec=(1+skip_first):(length(search_for)-skip_last);

a=jet;
col_L=length(search_for)-skip_first-skip_last;
if col_L==0;col_L=1;end
int=floor(length(a)/(col_L));
col_mat=a(1:int:int*(col_L),:);

for y=1:length(syl_vec)           % y indexes syllables, x indexes occurances

    y_val=syl_vec(y);

    x_data{y}=MAT_nspikes(:,y_val);xlab='# spikes, during syl';x_data_all=MAT_nspikes(:,syl_vec);

    y_data{y,1}=pitch(:,y_val);ylab{1}='pitch (Hz)';y_data_all=pitch(:,syl_vec);
    y_data{y,2}=log10(MAT_amplitude_16msec_around_pitchquant(:,y_val));ylab{2}='Mean LOG10 power 16msec around pitch quant (arb units)';y_data_all=MAT_amplitude_16msec_around_pitchquant(:,syl_vec);

    id_neginf=find(y_data{2}==-Inf);
    if ~isempty(id_neginf);disp('Setting -Inf amplitudes to zero');y_data{2}(id_neginf)=0;end

    y_data{y,3}=entropy(:,y_val);ylab{3}='entropy';
    y_data{y,4}=serial_order_tmp(:,y_val);ylab{4}='serial order';
    y_data{y,5}=harmonics(:,y_val);ylab{5}='hi/lo harmonics power ratio';
    y_data{y,6}=duration(:,y_val);ylab{6}='duration';
    y_data{y,7}=amp_slope(:,y_val);ylab{7}='amplitude slope';
%     y_data{y,8}=pitch_slope_tmp(:,y_val);ylab{8}='pitch slope';
%     y_data{y,9}=entropy_slope_tmp(:,y_val);ylab{9}='entropy slope';

    for z=1:length(y_data)
        mean_x{z}(y)=mean(x_data{y});
        mean_y{z}(y)=mean(y_data{z});

        if z<=4
            az(1)=subplot(3,4,z);hold on
            %        az(1)=subplot(3,3,z);hold on
            %        az(1)=subplot(3,2,z);hold on

            clear xlt
            xlt{1}='Number of spikes during specified interval';
            %xlt{1}=['# of spikes from [-' num2str(pre_onset) ' to ' num2str(post_onset) '] ms around onset']  ;
            %        xlt{2}=RemoveUnderScore(loadfilename);

            if spikes_on_x_axis
                plot(mean_x{z}(y),mean_y{z}(y),['k' symbol],'markerfacecolor',col_mat(y,:))
                xlabel(xlt)
                ylabel(ylab{z})
            else
                plot(mean_y{z}(y),mean_x{z}(y),['k' symbol],'markerfacecolor',col_mat(y,:))
                ylabel(xlt)
                xlabel(ylab{z})
            end
        end
    end

    %  plotting
    az(1)=subplot(3,4,1);
    %    az(1)=subplot(3,3,1);
    %    az(1)=subplot(3,2,1);

    mean_pitch_out(y)=mean_y{1}(y);
    mean_amp_out(y)=mean_y{2}(y);
    mean_entropy_out(y)=mean_y{3}(y);
    mean_serial_order_out(y)=mean_y{4}(y);
    mean_harmonics_out(y)=mean_y{5}(y);
    mean_durations_out(y)=mean_y{6}(y);
    mean_amp_slope_out(y)=mean_y{7}(y);
%     mean_pitch_slope_out(y)=mean_y{8}(y);
%     mean_entropy_slope_out(y)=mean_y{9}(y);
    mean_spikes_out(y)=mean_x{1}(y);
end

% structure array with mean output
mean_outs.pitch=mean_pitch_out;
mean_outs.amp=mean_amp_out;
mean_outs.entropy=mean_entropy_out;
mean_outs.harmonics=mean_harmonics_out;
mean_outs.spikes=mean_spikes_out;
mean_outs.serial_order=mean_serial_order_out;
mean_outs.duration=mean_durations_out;
mean_outs.amp_slope=mean_amp_slope_out;
% mean_outs.pitch_slope=mean_pitch_slope_out;
% mean_outs.entropy_slope=mean_entropy_slope_out;

function [x_data,y_data,regress_struct]=plot_and_regress(x_data,y_data,col_mat,col_rank,xlab,ylab)

% % NEW - remove outlier >x stdeviations greater than mean - this is for
% % occaisional ampltude outliers
if 0
    disp('NOT REMOVING OUTLIERS')
else    % remove outliers based on y data
    [x_data,y_data]=remove_outliers(x_data,y_data);
end

[FIT,tmp,tmp,tmp,STATS] =regress(y_data,[x_data ones(size(x_data))]); % fit 1 is slope, fit 2 is intercept
%Rsquared=.01*round(100*STATS(1));
Rsquared_str=num2str(STATS(1));
pval_str=num2str(STATS(3));
plot(x_data,y_data,'.','color',col_mat(col_rank,:));hold on
%plot(mean(x_data),mean(y_data),'ko','markerfacecolor',col_mat(col_rank,:));hold on
xl=get(gca,'xlim');
yline=FIT(1)*xl+FIT(2);
plot(xl,yline,'color',col_mat(col_rank,:));hold on
set(gca,'xlim',xl)
if ~isempty(xlab)
    xlabel(xlab)
end
if ~isempty(ylab)
    ylabel(ylab)
end
title(['R^2 = ' Rsquared_str ', p = ' pval_str ', slope = ' num2str(FIT(1))])
regress_struct.FIT=FIT;
regress_struct.STATS=STATS;


function [seq_array,skip_first,skip_last,t_assay]=generate_seq_array(birdname,syllable_to_quant)


if strcmp(birdname,'pu77bk41')
    %     if strcmp(syllable_to_quant,'a')
    %         skip_first=0;    skip_last=0;
    %         seq_array{1}='a';  % 'a'
    if strcmp(syllable_to_quant,'b')
        skip_first=0;    skip_last=1;
        seq_array{1}='bc';
    elseif strcmp(syllable_to_quant,'c')
        skip_first=0;    skip_last=0;
        seq_array{1}='c';
%     elseif strcmp(syllable_to_quant,'d')
%         skip_first=0;    skip_last=0;
%         seq_array{1}='d';
%     elseif strcmp(syllable_to_quant,'d')
%         skip_first=1;    skip_last=1;
%         seq_array{1}='cdd';
    elseif strcmp(syllable_to_quant,'d')
        skip_first=1;    skip_last=1;
        seq_array{1}='ddd';
    elseif strcmp(syllable_to_quant,'e')
        skip_first=0;    skip_last=1;
        seq_array{1}='ef';
    elseif strcmp(syllable_to_quant,'f')
        skip_first=1;    skip_last=0;
        seq_array{1}='ef';
    end
    disp('new calc of premotor window')
    [tmp,t_assay,tmp]=syllable_params_by_bird(birdname,syllable_to_quant(1));
elseif strcmp(birdname,'o85pu54')
    skip_first=0;    skip_last=0;
    seq_array{1}=syllable_to_quant;
    disp('new calc of premotor window');pause(1)
    [tmp,t_assay,tmp]=syllable_params_by_bird(birdname,syllable_to_quant(1));

elseif strcmp(birdname,'bl82bl81')
    if strcmp(syllable_to_quant,'d')
        skip_first=1;   skip_last=0; 
        seq_array{1}='dd';
    else
        skip_first=0;    skip_last=0;
        seq_array{1}=syllable_to_quant;
    end
    [tmp,t_assay,tmp]=syllable_params_by_bird(birdname,syllable_to_quant(1));
elseif strcmp(birdname,'r93bl81')
    if strcmp(syllable_to_quant,'a')
        skip_first=0;    skip_last=0;
        seq_array{1}='a';
    elseif strcmp(syllable_to_quant,'b')
        skip_first=0;    skip_last=0;
        seq_array{1}='b';
    elseif strcmp(syllable_to_quant,'c')
        skip_first=0;    skip_last=0;
        seq_array{1}='c';
    elseif strcmp(syllable_to_quant,'e')
        skip_first=1;    skip_last=1;
        seq_array{1}='eee';
    else
        skip_first=0;    skip_last=0;
        seq_array{1}=syllable_to_quant;
    end
    [tmp,t_assay,tmp]=syllable_params_by_bird(birdname,syllable_to_quant(1));
    % % elseif strcmp(birdname,'Pk35G27')
    % %     if strcmp(syllable_to_quant,'a')    % this discrimination added 12/6/2007
    % %         skip_first=1;    skip_last=1;
    % %         seq_array{1}='bad';
    % %     else
    % %     skip_first=0;    skip_last=0;
    % %         seq_array{1}=syllable_to_quant;
    % %     end
    % %     [tmp,t_assay,tmp]=syllable_params_by_bird_MEL(birdname,syllable_to_quant(1));

elseif strcmp(birdname,'B53O71') |  strcmp(birdname,'r93bl81') |  strcmp(birdname,'Pk35G27') |  strcmp(birdname,'o85pu54') |  strcmp(birdname,'G45G46') |  strcmp(birdname,'Pu55Pu22') |  strcmp(birdname,'W90O73') |  strcmp(birdname,'W32Pi51') |  strcmp(birdname,'W15W94')|  strcmp(birdname,'O14O15')|    strcmp(birdname,'W96Pi45')
    skip_first=0;    skip_last=0;
    seq_array{1}=syllable_to_quant;
    [tmp,t_assay,tmp]=syllable_params_by_bird(birdname,syllable_to_quant(1));
elseif strcmp(birdname,'g38o18') | strcmp(birdname,'g38o18')   | strcmp(birdname,'p85g54')  | strcmp(birdname,'g91pu54') 
    skip_first=0;    skip_last=0;
    seq_array{1}=syllable_to_quant;
    [tmp,t_assay,tmp]=syllable_params_by_bird(birdname,syllable_to_quant(1));
elseif strcmp(birdname,'b39b14') | strcmp(birdname,'b9r63')  | strcmp(birdname,'blueblue')  | strcmp(birdname,'jman3')  | strcmp(birdname,'masa') | strcmp(birdname,'pur98w3') % mimi lman ZF
    skip_first=0;    skip_last=0;
    seq_array{1}=syllable_to_quant;
    [tmp,t_assay,tmp]=syllable_params_by_bird(birdname,syllable_to_quant(1));
else
    error('Unknown birdname')
end
if ~exist('t_assay');error('Unknown syllable');end


% remove outliers based on y data
function [x_data,y_data]=remove_outliers(x_data,y_data);

% disp('SKIPPING OUTLIER REMOVAL STEP')
% return

L_orig=length(y_data);mean_y=mean(y_data);std_y=std(y_data);
n_std=4;
cutoff_hi=mean_y+n_std*std_y;
cutoff_lo=mean_y-n_std*std_y;

id=find(y_data<=cutoff_hi & y_data>=cutoff_lo);
if length(id)<L_orig
    %    disp('  ')
    disp(['Deleting ' num2str(L_orig-length(id)) ' outliers > ' num2str(n_std) ' stds beyond mean'])
    %    disp('  ')
end
y_data=y_data(id);
x_data=x_data(id);

set(gcf,'paperposition',[0.3562    5.9531    4.3875    4.9906])