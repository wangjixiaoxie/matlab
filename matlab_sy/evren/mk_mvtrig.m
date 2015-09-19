function mk_mvtrig(batch);
% usage : mk_mvtrig(batch);
%
% takes in a file with a list of filenames, and makes a 
% shell script to move those files 
% this is necessary because moving them from matlab takes forever
% and i have not had a chance to write a perl script to do it yet

fid = fopen(batch,'r');
fid2 = fopen('mvtrig','w');
fprintf(fid2,'#!/bin/sh\n');
fprintf(fid2,'mkdir triggers\n');
while(~feof(fid))
	fn = fscanf(fid,'%s',1);
	
	if length(fn>2)
		p = findstr(fn,'.');
		if (length(p)==1)
			cmd = ['mv ',fn(1:p(end)),'* triggers/'];
			fprintf(fid2,'%s\n',cmd);
		else
			cmd = ['mv ',fn(1:p(end)),'* triggers/'];
			fprintf(fid2,'%s\n',cmd);
		end
	end

end
fclose(fid);fclose(fid2);
eval(['!chmod u+x mvtrig']);
return;
