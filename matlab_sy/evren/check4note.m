function cnt=check4note(bt,NT);
%

cnt=0;
fid=fopen(bt,'r');
while (1)
	fn=fgetl(fid);
	if (~ischar(fn));
		break;
	end
	if (~exist([fn,'.not.mat']))
		disp(['File ',fn,'.not.mat does not exist']);
		continue;
	end
	load([fn,'.not.mat']);
	if (length(findstr(labels,NT))==0)
		disp(fn);
		cnt=cnt+1;
	end
end
fclose(fid);
return;
