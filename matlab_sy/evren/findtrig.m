function [] = findtrig(batchfile,filetype,thresh);
% findtrig(batchfile,filetype,thresh);
%
% goes through the batch file of obs1r files looking for 
% BirdTAF triggers
%

ftrigger = fopen([batchfile,'.trigs'],'a');
if (~exist('thresh'))
	thresh = 0.50;
end

FIRSTRUN=1;
fid = fopen(batchfile,'r');
while (1)
	fn  = fgetl(fid);
	if (~ischar(fn))
		break;
	end
	disp(fn);
	if (exist(fn,'file'))
		[taf_out,fs]=evsoundin('',fn,filetype);
	end
	
	if (length(taf_out)>0)
		p = find(taf_out>thresh);
	
		if (length(p) > 0)
			if (FIRSTRUN==1)
				fprintf(ftrigger,'%s',fn);
				FIRSTRUN=0;
                        else
				fprintf(ftrigger,'\n%s',fn);
                        end
	         end
	end
end
fclose(fid);fclose(ftrigger);
return;
