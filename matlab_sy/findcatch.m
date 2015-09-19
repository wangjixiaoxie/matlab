function findcatch(batch);
% findcatch(batch);
% write all file names in batch file that are catch trials
% into [batch,'.catch'];
fid=fopen(batch,'r');
fcat=fopen([batch,'.catch'],'w');
fnot=fopen([batch,'.notcatch'],'w');
while (1)
	fn=fgetl(fid);
	if (~ischar(fn))
		break;
	end
	if (~exist(fn,'file'))
		continue;
	end
	pp=findstr(fn,'.');
	if (length(pp)<1)
		continue;
	end
	fnr = [fn(1:pp(end)),'rec'];
	if (~exist(fnr,'file'))
		disp(['Rec file does not exist = ',fnr]);
		continue;
	end
	frec=fopen(fnr,'r');
	while (1)
		ln = fgetl(frec);
		if (~ischar(ln))
			break;
		end
		pp = findstr(ln,'Catch');
		if (length(pp)==1)
			pp = findstr(ln,'=');
			iscatch = str2num(ln((pp(end)+1):end));
			if (iscatch==1)
				fprintf(fcat,'%s\n',fn);
            else
                fprintf(fnot,'%s\n',fn);
			end
		end
	end
	fclose(frec);
end
fclose(fid);fclose(fcat);fclose(fnot);
return;
