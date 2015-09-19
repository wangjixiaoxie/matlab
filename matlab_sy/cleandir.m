function cleandir(batchnotmat);
%cleandir(batchnotmat);
%
fid  = fopen(batchnotmat,'r');
fout = fopen('rmdata','w');
fprintf(fout,'#!/bin/sh\n');
while (1)
	fn = fgetl(fid);
	if (~ischar(fn))
		break;
	end
	if (~exist(fn,'file'))
		continue;
	end
	load(fn);

	if (length(findstr(labels,'i'))==length(labels))
		pp = findstr(fn,'.not.mat');
		fprintf(fout,'rm %s*\n',fn(1:(pp-1)));
	end
end
fclose(fid);fclose(fout);
eval(['!chmod u+x rmdata']);
return;
