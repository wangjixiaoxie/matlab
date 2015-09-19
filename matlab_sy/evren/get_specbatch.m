function get_specbatch(bt,outfile);
%get_specbatch(bt,outfile);

fid=fopen(bt,'r');
if (fid==-1)
    disp(['Can''t find batch file : ',bt]);
    return;
end

fid2=fopen([bt,'.out'],'w');
while (1)
        
    fn=fgetl(fid);
    if (~ischar(fn))
        break;
    end
    if (~exist(fn,'file'))
        continue;
    end
    
    rd=readrecf(fn);
    if (strcmp(rd.outfile,outfile))
        fprintf(fid2,'%s\n',fn);
    end
end
fclose([fid,fid2]);
return;