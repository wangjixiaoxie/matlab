function vals=get_reltime(bt,REFNT,NT,PRENT,PSTNT);
% vals=get_reltime(bt,REFNT,NT,PRENT,PSTNT);
%

if (~exist('PRENT','var'))
	PRENT='';
elseif (length(PRENT)==0)
	PRENT='';
end

if (~exist('PSTNT','var'))
	PSTNT='';
elseif (length(PSTNT)==0)
	PSTNT='';
end

if (length(findstr([PRENT,NT,PSTNT],REFNT))==0)
	disp(['REFNT must be in search str']);
	disp(['Search string = ''',PRENT,NT,PSTNT,'''']);
	disp(['Ref note      = ''',REFNT,'''']);
	return;
end
refind=findstr([PRENT,NT,PSTNT],REFNT);
refind=refind(1);

vals=[];
ff=load_batchf(bt);
for ii=1:length(ff)
	fn=[ff(ii).name,'.not.mat'];
	if (~exist(fn,'file'))
		continue;
	end
	load(fn);

	pp=findstr(labels,[PRENT,NT,PSTNT])+length(PRENT);
	if (length(pp)==0)
		continue;
	end
	for jj=1:length(pp)
		refnt_t=onsets(pp(jj)-length(PRENT)+refind-1);
		nt_t=onsets(pp(jj));
        	vals=[vals;nt_t-refnt_t];
	end
end
return;
