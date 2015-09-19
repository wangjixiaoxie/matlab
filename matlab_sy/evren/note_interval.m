function [vals,intervals]=note_interval(batchnotmat,pattern,ADDNOTMAT);
% [vals,intervals]=note_interval(batchnotmat,pattern,ADDNOTMAT);
%
if (~exist('ADDNOTMAT'))
	ADDNOTMAT = 0;
end

vals = [];
fid = fopen(batchnotmat);cnt=0;
while (1)
	fn = fgetl(fid);
	if (~ischar(fn))
		break;
	end
	if (ADDNOTMAT)
		fn = [fn,'.not.mat'];
	end
	if (~exist(fn,'file'))
		continue;
	end
	load(fn);

	cnt = cnt + 1;
	intervals(cnt).fn = fn;
	intervals(cnt).pat= pattern;
	pos = findstr(labels,pattern);
	if (length(pos)>0)
		tmpint = zeros([length(pos),length(pattern)-1]);
		for ii = 1:length(pos)
			for jj = 1:length(pattern)-1
				tmpint(ii,jj) = ...
				onsets(pos(ii)+jj)-offsets(pos(ii)+jj-1);
				%onsets(pos(ii)+jj)-onsets(pos(ii)+jj-1);
			end
		end
		intervals(cnt).int = tmpint;
		vals = [vals;tmpint];
	else
		intervals(cnt).int = [];
	end
end
fclose(fid);

return;
