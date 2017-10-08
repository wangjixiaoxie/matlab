%% to change name of classes - i.e. add numbers, so can combine with other analysis (e./g. if crashed)
% e.g. go to Results_xaaa_AlgnSyl2Onset1_27Sep2017_0245_RAallbirds
fnames = dir('classes*');

for i=1:length(fnames)
   
    dotloc = strfind(fnames(i).name, '.mat');
    
    iternum = fnames(i).name(8:dotloc-1);
    
    eval(['!mv classes' num2str(iternum) '.mat classes88' num2str(iternum) '.mat']);
    eval(['!mv params' num2str(iternum) '.mat params88' num2str(iternum) '.mat']);
    
    
end