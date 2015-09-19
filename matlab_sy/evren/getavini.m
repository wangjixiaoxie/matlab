function vals=getavini(batch,NT);
%vals=getavini(batch,NT);

vals=[];
fid=fopen(batch,'r');
while (1)
	fn=fgetl(fid);
	if (~ischar(fn));
		break;
	end
	if (~exist([fn,'.not.mat']))
		continue;
	end

	load ([fn,'.not.mat']);
	pp=findstr(labels,NT);
	if (length(pp)<2)
		continue;
	end
	for ii = 1:length(pp)-1
		dt = onsets(pp(ii+1))-onsets(pp(ii));
		vals=[vals;dt];
	end
end
fclose(fid);
return;
