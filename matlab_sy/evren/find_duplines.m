function find_duplines(filename);
% find_duplines(filename);

fid=fopen(filename,'r');
fid2=fopen([filename,'.nodup'],'w');
ln2=fgetl(fid);
fprintf(fid2,'%s\n',ln2);
ln1=fgetl(fid);
while (1)
	if (~ischar(ln1))
		break;
	end
	if (strcmp(ln1,ln2))
		disp(ln1);
    else
        fprintf(fid2,'%s\n',ln1);
	end
	ln2=ln1;
	ln1=fgetl(fid);
end
fclose([fid,fid2]);
return;
	
