function new=jwresamp(old,oldrate,newrate);
% usage : new=jwresamp(old,oldrate,newrate)
%if every 5 indices of old vector will become 2 of the new, oldrate=5,newrate=2

fid = fopen(batch,'r');
fid2 = fopen('rmdata','w');
fprintf(fid2,'#!/bin/sh\n');
while(1)
	fn = fgetl(fid);
    if (~ischar(fn))
        break;
    end
	
	if length(fn>2)
		p = findstr(fn,'.');
		if (length(p)>0)
            if (onlyebin==1)
                cmd = ['rm ',fn];
                fprintf(fid2,'%s\n',cmd);
            else
                cmd = ['rm ',fn(1:p(end)),'*'];
                fprintf(fid2,'%s\n',cmd);
            end
		else
			cmd=['rm ',fn];
			fprintf(fid2,'%s\n',cmd);
		end
	end

end
fclose(fid);fclose(fid2);
eval(['!chmod u+x rmdata']);
return;
