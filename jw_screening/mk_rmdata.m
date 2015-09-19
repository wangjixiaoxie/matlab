function mk_rmdata(batch,onlyebin);
% usage : mk_rmdata(batch,rec_too);
%% EXAMPLE mk_rmdata('bk49bk25.del',1);
% then go into window and type ./rmdata
% set onlybin=0 to also delete the .rec files
%
% takes in a file with a liost of filenames, and makes a 
% shell script to delete those files
% this is necessary because deleting them from matlab takes forever
% and i have not had a chance to write a perl script to do it 


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
