function [score, minsfound]=testing(smoothed)
%[score, minsfound]=mpeaks(smoothed)
%score is peak-multiplicity score. minsfound is a list of indices (in samps) into the smoothed syl vector that are local minima.


slopes=diff(smoothed);
downs=slopes<0;
minmax=conv([1 -1],downs);
find(minmax) %only for debugg-- compare with after deleting close mms

if minmax(length(minmax)) ~= -1
	disp('Smoothed syl didn''t end in a local min.. calling the end a minimum.')
	minnmax(length(minmax)) = -1;
end
minmax(1)=-1;

%to delete close pairs:
mm_inds=find(minmax);
for k=1:length(mm_inds)
smoothed_mms(k)=smoothed(mm_inds(k));
end
closemms=find(abs(diff(smoothed_mms))<1000); %that's indices into smoothed_mms, which has each element corresponding through mm_inds to a point in minmax (or smoothed).
if length(closemms)>=1
for k=1:length(closemms)
	minmax(mm_inds(closemms(k))) %to check if I get the right #
	minmax(mm_inds(closemms(k)))=0
	minmax(mm_inds(closemms(k+1)))=0
	mm_inds=find(minmax);
	for p=1:length(mm_inds)
		smoothed_mms(p)=smoothed(mm_inds(p));
	end
	
end
end
find(minmax) %compare to before 'cleanup'
