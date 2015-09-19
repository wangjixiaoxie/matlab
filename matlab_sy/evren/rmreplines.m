function rmreplines(infile,outfile);

fid=fopen(infile,'r');
fid2=fopen(outfile,'w');
fn=fgetl(fid);
fprintf(fid2,'%s',fn);
cnt=0;
while (1)
	fn2=fgetl(fid);
	if (~ischar(fn2))
		break;
	end
	if (~strcmp(fn2,fn))
		fprintf(fid2,'\n%s',fn2);
	else	
		cnt=cnt+1;
	end
	fn=fn2;
end
fclose([fid,fid2]);
disp(['N=',num2str(cnt)]);
return;
