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
end
fclose(fid);fclose(fcat);fclose(fnot);
return;
