function [peakpos,pkheight,scores,lengths,labs]=colorscale(n,peakpos,pkheight,scores,lengths,labs)
%[peakpos,pkheight,scores,lengths,labs]=colorscale(n, peakpos,pkheight,scores,lengths,labs);
%adds a final element to 5 vectors:
%labs, pkheight, peakpos, lengths, and scores.
%for all but labs, final element is the max of that vector.
%for labs it is n.
%checks if it looks like colorscale has already been used before,
%if so it does not add elements to the vectors, it just resets the
%last element.


if (peakpos(length(peakpos))==max(peakpos))&(lengths(length(lenghts))==max(lengths))&(pkheight(length(pkheight))==max(pkheight))&...
(scores(length(scores))==max(scores))&(0==any(find(labs(1:length(labs)-1)==labs(length(labs)))))
	disp('colorscale seems to have been used on these before.. redoing.')
	i=0;
	else
	disp('all five vectors will get a new element appended.')
	i=1;
end
peakpos(length(peakpos)+i)=max(peakpos);
lengths(length(lengths)+i)=max(lengths);
pkheight(length(pkheight)+i)=max(pkheight);
scores(length(scores)+i)=max(scores);
labs(length(labs)+i)=n;
