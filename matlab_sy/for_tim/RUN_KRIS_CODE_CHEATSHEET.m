function foo
!dir /B *shift*.cbin > batchfile

autolabel_KB('batchfile','obs2r','foo')

load foo_tmpstrct.template_strct.mat
get_mahal_syl(foo_tmpstrct,'batchfile','foo2',1)


if 0
% run evsonganaly
% set thresh low
% will show 10 examples of each syl, starting with center exemplar.  if
% first one is good example, pick it, and ONLY label that one (one of each
% syl).

make_avg_spects('batchfile','m','foo2.mhlsyl.min_strct.mat.mahal_syls.cbin','obs0','foo3') % just changed from obs0r to obs0

%load foo3FLOOR.spct_strct.mat
load foo3.spct_strct.mat
plot_spct_strct_SAM_HACK(spct_strct)

plot_spct_strct(spct_strct)
end

if 0
    clear
    dir_arr{1}='E:\pitch_shift\r12r11\03_30_2007';
    dir_arr{2}='E:\pitch_shift\r12r11\05_29_2007';
%    dir_arr{3}='E:\pitch_shift\r12r11\03_16_2007';
    for x=1:length(dir_arr)
        cd(dir_arr{x})
        disp(pwd)
        RUN_KRIS_CODE_CHEATSHEET
    end
    
end