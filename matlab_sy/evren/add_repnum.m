for ii=1:length(fv)
	labels=fv(ii).lbl;
	ind=fv(ii).ind;
	NT=labels(ind);
	kk=ind-1;cnt=1;
	while (kk>0)
		if (~strcmp(NT,labels(kk)))
			break;
		else
			kk=kk-1;
			cnt=cnt+1;
		end
	end
	fv(ii).repnum=cnt;
	
	kk=ind;
	while (1)
		if (kk>=length(labels))
			break;
		end
		if (strcmp(labels(kk+1),NT))
			kk=kk+1;
		else
			break;
		end
	end
	fv(ii).reptot=(kk-ind)+cnt;
end
