function vals=get_durs(batch,NT,PRENT,PSTNT);
% vals=get_durs(batch,NT,PRENT,PSTNT);
% vals=[datenum,dur,filenum,noteindex];


ff=load_batchf(batch);
vals=[];TH=0;str=[PRENT,NT,PSTNT];
for ii=1:length(ff)
	if (~exist([ff(ii).name,'.not.mat']))
		continue;
	end

	load([ff(ii).name,'.not.mat']);
	pp=findstr(labels,'0');labels(pp)='-';
	if (TH==0)
		TH=threshold;
	else
		if (TH~=threshold)
			disp(['Hey ',fn])
			TH=threshold;
		end
	end

	dn=fn2datenum(ff(ii).name);
	pp=findstr(labels,str).';pp=pp+length(PRENT);
	vals=[vals;dn*ones(size(pp)),offsets(pp)-onsets(pp),ii*ones(size(pp)),pp];
end
return;
