function [smvals,refntvals]=smanaly(bt,NTS,REFNT,CS);
% smvals=smanaly(bt,NTS,REFNT,CS);
%
smvals=[];
for kk=1:length(NTS)
	smvals(kk).vals=[];
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

	sm=evsmooth(dat,fs,0.0);
	sm=sqrt(sm);

	inds=[floor(onsets(1)*1e-3*fs),ceil(offsets(end)*1e-3*fs)];
	normv=mean(sm([inds(1):inds(2)]));
	sm=sm./normv;

	for kk=1:length(NTS)
		NT=NTS(kk).NT;
		PRENT=NTS(kk).PRENT;
		PSTNT=NTS(kk).PSTNT;
		TOFF =NTS(kk).TOFF;
		pp=findstr(labels,[PRENT,NT,PSTNT])+length(PRENT);
		if (length(pp)==0)
			continue;
		end

		for jj=1:length(pp)
			inds=[0,0];
			if (TOFF==0)
			    inds=[floor(onsets(pp(jj))*1e-3*fs),...
		      	          ceil(offsets(pp(jj))*1e-3*fs)];
			else
			    inds(1)=floor((onsets(pp(jj))*1e-3+TOFF)*fs);
			    inds(2)=inds(1)+ceil(NTS(kk).DT*fs);
			end
			newval=sum(sm([inds(1):inds(2)]))./fs;
		
			smvals(kk).vals=[smvals(kk).vals,newval];
		end
	end
end

for kk=1:length(NTS)
	smvals(kk).mn=mean(smvals(kk).vals);
	smvals(kk).st=std(smvals(kk).vals);
end
smvals=rmfield(smvals,'vals');

return;
