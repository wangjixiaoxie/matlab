function label_othernotes(bt,NT,times,newnts);
% label_othernotes(bt,NT,times,newnts);
% bt batch files name (.cbin files)
% NT is the already labeled NT
% times are in ms
% newnts are the notes to put in
%

if (length(times)~=length(newnts))
	disp(['ERROR: times and newnts do not have the same size']);
	return;
end

ff=load_batchf(bt);

for ii=1:length(ff)
	fn=ff(ii).name;
	if (exist([fn,'.not.mat']))
		load([fn,'.not.mat']);
		disp(fn);
	else
		disp(['no file : ',fn,'.not.mat']);
		continue;
	end

	pp=findstr(labels,NT);

	if (length(pp)==0)
		continue;
	end

	for jj=1:length(pp)
		ind=pp(jj);
		for kk=1:length(times)
			lbl=newnts(kk);
			tt=times(kk)+onsets(pp(jj));
			ppp=find((onsets<=tt)&(offsets>=tt));
			if (length(ppp)>0)
				if (ppp>ind)
					labels(ppp)=lbl;
				end
			end
		end
	end
	eval(['save -append ',fn,'.not.mat labels']);
end

return;
