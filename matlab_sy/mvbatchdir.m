function mvbatchdir(batch,dirname)
    
    fid=fopen(batch,'r');
    
while (1)
	fn=fgetl(fid);
	if (~ischar(fn))
		break;
	end
	if (~exist(fn,'file'))
		continue;
	end
	disp(fn)
	
    strcmd=['!mv ' fn ' ' dirname];
    eval(strcmd)
end
fclose(fid);
return