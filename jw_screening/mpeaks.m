function [score, minsfound, maxsfound]=mpeaks(smoothed)
%[score, minsfound, maxsfound]=mpeaks(smoothed)
%score is peak-multiplicity score. minsfound is a list of indices (in samps) into the smoothed syl vector that are local minima.


slopes=diff(smoothed);
downs=slopes<0;
minmax=conv([1 -1],downs);
%%%%%%%%%%%%%%%%%%%%

%find(minmax) %only for debugg-- compare with after deleting close mms

if minmax(length(minmax)) ~= -1
	disp('Smoothed syl didn''t end in a local min.. calling the end a minimum.')
	minnmax(length(minmax)) = -1;
end
minmax(1)=-1;
%%%%%%%%test, remove this segment%%%%%%

smoothed(1)=smoothed(1)/2;
smoothed(length(smoothed))=smoothed(length(smoothed))/2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%to delete close pairs:
%mm_inds=find(minmax);
%for k=1:length(mm_inds)
%smoothed_mms(k)=smoothed(mm_inds(k));
%end
%closemms=find(abs(diff(smoothed_mms))<1000); %that's indices into smoothed_mms, which has each element corresponding through mm_inds to a point in minmax (or smoothed).
%if length(closemms)>=1
%for k=1:length(closemms)
%	minmax(mm_inds(closemms(k))) %to check if I get the right #
%	minmax(mm_inds(closemms(k)))=0
%	minmax(mm_inds(closemms(k+1)))=0
%	mm_inds=find(minmax);
%	for p=1:length(mm_inds)
%		smoothed_mms(p)=smoothed(mm_inds(p));
%	end
	
%end
%end
%find(minmax) %compare to before 'cleanup'


%minmax(1)=-1;

%%%%%%%%%%%%%%%%%%%%%

mins=find(minmax<0);
maxs=find(minmax>0);
d=0;
score=0;
for i=1:(length(mins)-1)		%these are initials and
	for f=0:(length(mins)-i-1)	%finals -- the beginning and end values that define this domain.
	%%%%%%%%%%%%testing--delete this segment%%%%%%%%
	
	%skip any domain that is not "significant"
%	if max(smoothed(mins(i):mins(length(mins)-f)))-smoothed(mins(i))<.05*mean(smoothed)|max(smoothed(mins(i):mins(length(mins)-f)))-smoothed(mins(length(mins)-f))<.05*mean(smoothed)
%	disp(['skipping domain', num2str(mins(i)), ' to ',num2str(mins(length(mins)-f))])
%	break
%	end
	
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		d=d+1;
		domains(d,:)=[mins(i), mins(length(mins)-f), 666, 666, 666, 666];
		%first, find this domain's score increment
		imax=find(maxs>domains(d,1));
		imax=maxs(imax(1));	%had to do that in 2 steps because find(maxs>domains(d,1), 1, 'first') wouldn't work
		fmax=find(maxs<domains(d,2));
		fmax=maxs(fmax(length(fmax)));
		domincr=(min((smoothed(imax)-smoothed(domains(d,1))),(smoothed(fmax)-smoothed(domains(d,2))))/max((smoothed(imax)-smoothed(domains(d,1))),(smoothed(fmax)-smoothed(domains(d,2)))));
		%if min((smoothed(imax)-smoothed(domains(d,1))),(smoothed(fmax)-smoothed(domains(d,2))))==(smoothed(imax)-smoothed(domains(d,1)))
		%smaller='i'
		%else
		%smaller='f'
		%end
		%then, subtract this domain's decrement
		valvals=[];
		if (i+1)<(length(mins)-f)
			for v=1:length(mins)-f-i-1
				valvals(v)=smoothed(mins(i+v));
			end
			domdecr=(min(smoothed(imax), smoothed(fmax))-min(valvals))/(min(smoothed(imax), smoothed(fmax)));

		else
		domdecr=0;
		%no decrement of there are no valleys w/in the domain
		end
		%if domscore<0
		%domscore=0;
		%end
		domscore=domincr-(domincr*domdecr);
		%%bullshit factor below!!
		%%if domscore<0
		%%domscore=.2*domscore;
		%%end
		score=score+domscore;
		domains(d,3)=domincr;domains(d,4)=domdecr;
		%disp([d, domains(d,:), valvals])
	end
end

minsfound=mins;
maxsfound=maxs;
%domains(1:10,:)

