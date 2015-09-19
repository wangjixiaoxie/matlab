function [nlen]=evnotestat(batchnotmat,note);
%
nlen=[];
fid=fopen(batchnotmat,'r');
while (1)
	fn=fgetl(fid);
	if (~ischar(fn))
		break;
	end

	if (~exist(fn,'file'))
		break;
	end
	load(fn);

	p = findstr(lower(labels),note);
	if (length(p)==0)
		continue;
	end
	if (length(p)>1)
		nlennew=[];
		ind1 = 1;ind2=1;
		while (ind1<=length(p))
			while (ind2<length(p))
				if (p(ind2+1)==(p(ind2)+1))
					ind2 = ind2 +1;
                                else
					break;
			        end
		        end
			nlennew = [nlennew;offsets(p(ind2))-onsets(p(ind1))];
			ind1=ind2+1;ind2=ind1;
		end
	else
		nlennew = [offsets(p)-onsets(p)];
	end
	nlen=[nlen;nlennew];
end
fclose(fid);
return;
