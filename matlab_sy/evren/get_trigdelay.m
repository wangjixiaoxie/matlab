function vals=get_trigdelay(batch,NT);
% vals=get_trigdelay(batch,NT);
%

if (~exist('NT'))
	NT=[];
end

vals=[];
fid=fopen(batch,'r');
while (1)
	fn=fgetl(fid);
	if (~ischar(fn));
		break;
	end

	if (~exist([fn,'.not.mat'],'file'))
		continue;
	end
	load([fn,'.not.mat']);
	rd=readrecf(fn);

	tt=rd.ttimes;
	if (length(tt)==0)
		continue;
	end

	for ii=1:length(tt)
		pp=find((onsets<=tt(ii))&(offsets>=tt(ii)));
		if (length(pp)>0)
			if (length(NT)>0)
				if (strcmp(NT,labels(pp(1))))
					vals=[vals;...
                                   	 tt(ii)-onsets(pp(1)),...
                                         offsets(pp(1))-onsets(pp(1)),...
                                         fix(labels(pp(1)))];
				end
			else
				vals=[vals;...
                                      	tt(ii)-onsets(pp(1)),...
					offsets(pp(1))-onsets(pp(1)),...
                                        fix(labels(pp(1)))];
			end
		end
	end
end
fclose(fid);
return;
