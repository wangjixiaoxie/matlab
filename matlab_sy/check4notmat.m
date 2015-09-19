function check4notmat(batch,NOTE,DOLOWERCASE)
%check4notmat(batch,NOTE)

if (~exist('NOTE'))
	% only searchf or the existance of a .not.mat file
	NOTE = '';
end

if (~exist('DOLOWERCASE'))
	DOLOWERCASE=0;
end

fid =fopen(batch,'r');
fid2=fopen([batch,'.nonot'],'w');
while (1)
	fn=fgetl(fid);
	if (~ischar(fn))
		break;
	end
	pp=findstr(fn,'.not.mat');
	if (length(pp)>0)
		disp(['you''re batch file has .not.mat files?']);
		disp(['removing the .not.mat']);
		fn = fn(1:pp(end)-1);
	end
	disp(fn);

	if (~exist([fn,'.not.mat'],'file'))
		fprintf(fid2,'%s\n',fn);
	else
		load([fn,'.not.mat']);
		if (DOLOWERCASE==1)
			labels = lower(labels);
		end
		pp=findstr(labels,NOTE);
		if (length(pp)<1)
			fprintf(fid2,'%s\n',fn);
		end
	end
end
fclose(fid);fclose(fid2);
return;
