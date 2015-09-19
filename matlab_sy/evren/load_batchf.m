function ff=load_batchf(bt);
%returns a struct with batch file names
ff=[];
fid=fopen(bt,'r');
while (1)
	fn=fgetl(fid);
	if (~ischar(fn));
		break;
	end
	ff(length(ff)+1).name=fn;
end
fclose(fid);
return;
