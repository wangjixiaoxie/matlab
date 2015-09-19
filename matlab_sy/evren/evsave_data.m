function evsave_data(bfile);
% moves the data to keep into a temp folter
% then removes all .cbin in main directoy
% movbes files back after that
% used to purge out screened files which were not songs

input(['Data is going to be ERASED - proceed? press return']);

fid=fopen(bfile,'r');
cmd = ['!mkdir tmp_save_data'];eval(cmd);
while (~feof(fid))
	fn = fscanf(fid,'%s',1);
	p = findstr(fn,'.cbin');
	if (length(p)==1)
		cmd = ['!mv ',fn(1:p(1)),'* tmp_save_data/'];eval(cmd);
	end
end
fclose(fid);
eval(['!rm *.cbin']);
eval(['!rm *.rec']);
eval(['!mv tmp_save_data/* .']);
eval(['!rmdir tmp_save_data']);
return;
