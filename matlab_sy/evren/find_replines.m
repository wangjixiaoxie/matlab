function linenum=find_replines(bt);
%

fid = fopen(bt,'r');
fn=fgetl(fid);cnt=1;tmp=[];
while (1)
    if (~ischar(fn))
        break;
    end
    fn2=fgetl(fid);
    if strcmp(fn,fn2)
        tmp=[tmp;cnt];
        disp(fn);
    end
    cnt=cnt+1;
    fn=fn2;
end
linenum=tmp;
fclose(fid);
return; 
 