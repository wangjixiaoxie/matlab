function purge_data(bfile,rm_recfile);
% removes data in file bfile
% used to purge out screened files which were not songs

input(['Data is going to be ERASED - proceed? press return']);


fid=fopen(bfile,'r');
while (~feof(fid))
	fn = fscanf(fid,'%s',1)
	if (rm_recfile ~= 1)
		cmd = ['!rm ',fn]%;eval(cmd);
	else
		cmd = ['!rm ',fn];eval(cmd);
		p = findstr(fn,'.cbin');
		if (length(p)==1)
			cmd = ['!rm ',fn(1:p(1)),'rec'];eval(cmd);
		end
	end
end
fclose(fid);
return;
