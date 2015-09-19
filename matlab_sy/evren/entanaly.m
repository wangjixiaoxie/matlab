function entvals=entanaly(bt,NTS,REFNT,CS);
% entvals=smanaly(bt,NTS,REFNT,CS);
%
entvals=[];
for kk=1:length(NTS)
	entvals(kk).vals=[];
end

ff=load_batchf(bt);
for ii=1:length(ff)
	fn=ff(ii).name;
	if (~exist(fn,'file'))
		continue;
	end

	disp(fn);
	[dat,fs]=evsoundin('',fn,CS);
	load([fn,'.not.mat']);

	[sm,sp,t,f]=evsmooth(dat,fs,0.0);
	ppp=find((f>=500)&(f<=1e4));
	sp=abs(sp(ppp,:));
	sp=sp./(ones([size(sp,1),1])*sum(sp));
	ent=sum(-sp.*log2(sp));

	for kk=1:length(NTS)
		NT=NTS(kk).NT;
		PRENT=NTS(kk).PRENT;
		PSTNT=NTS(kk).PSTNT;
		TOFF =NTS(kk).TOFF;
		DT   =NTS(kk).DT;
		pp=findstr(labels,[PRENT,NT,PSTNT])+length(PRENT);
		if (length(pp)==0)
			continue;
		end

		for jj=1:length(pp)
			inds=[0,0];
			if (TOFF==0)
				[y,i]=min(abs(onsets(pp(jj))*1e-3-t));
				inds(1)=i;
				[y,i]=min(abs(offsets(pp(jj))*1e-3-t));
				inds(2)=i;
			else
				[y,i]=min(abs(onsets(pp(jj))*1e-3+TOFF-t));
				inds(1)=i;
				[y,i]=min(abs(offsets(pp(jj))*1e-3+TOFF+DT-t));
				inds(2)=i;
			end
			newval=mean(ent(inds(1):inds(2)));
		
			entvals(kk).vals=[entvals(kk).vals,newval];
		end
	end
end

for kk=1:length(NTS)
	entvals(kk).mn=mean(entvals(kk).vals);
	entvals(kk).st=std(entvals(kk).vals);
end
entvals=rmfield(entvals,'vals');

return;
