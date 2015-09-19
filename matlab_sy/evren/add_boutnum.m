for ii=1:length(fv)
	labels=fv(ii).lbl;
	ind=fv(ii).ind;
	NT=labels(ind);

	[tmp,vals]=get_repeats(labels);
	pp=find((vals(:,3)==fix('-'))|(vals(:,3)==fix('0')));
	vals(pp,:)=[];

	pp=find((ind>=vals(:,1))&(ind<=vals(:,2)));
	if (length(pp)>0)
		fv(ii).boutnum=pp(1);
	else
		fv(ii).boutnum=-1;
	end

	pp=findstr(labels,'-');
	labels(pp)=[];
	pp=findstr(labels,'0');
	labels(pp)=[];

	fv(ii).startbout=labels(1);
end
