function [outcell] = lt_batchsong_NamesToCell(batchf);
%% lt 8/17/16 - input batchfile, output cell with each line as entry

outcell={};
fid=fopen(batchf);
fline=fgetl(fid);

while ischar(fline)
    
    outcell=[outcell fline];
    
    
   fline=fgetl(fid);
end    