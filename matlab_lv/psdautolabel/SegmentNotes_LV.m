function [onsets, offsets]=SegmentNotes_LV(smooth, Fs, min_int, min_dur, threshold,lowerthreshold)
% [ons,offs]=evsegment(smooth,Fs,min_int,min_dur,threshold);
% segment takes smoothed filtered song and returns vectors of note
% onsets and offsets values are in seconds

h=[1 -1];

%threshold input
% notetimes=abs(diff(smooth))>threshold;
notetimes=smooth>threshold;
notetimes=single(notetimes);

%extract index values for note onsets and offsets
trans=conv(h,notetimes);
t_onsets  = find(trans>0);
t_offsets = find(trans<0);

onsets = t_onsets;offsets=t_offsets;
if ((length(onsets)<1)|(length(offsets)<1))
	onsets=[];offsets=[];
	return;
end

%same thing for lower threshold
lowernotetimes=smooth>lowerthreshold;
lowernotetimes=single(lowernotetimes);


if (length(t_onsets) ~= length(t_offsets))
    disp('number of note onsets and offsets do not match')
else
    
    %extract index values for note onsets and offsets
    lowertrans=conv(h,lowernotetimes);
    tempstart = 1;
%     tempend = t_onsets(2);
    for i = 1:length(t_onsets)
        if i<length(t_onsets)
            tempend = t_onsets(i+1);
        else
            tempend = length(lowertrans);
        end
        findlowon = find(lowertrans(tempstart:t_onsets(i))>0,1,'last');
        if ~isempty(findlowon)
         lowonset(i) = findlowon+tempstart;
        else
            lowonset(i) = nan;
        end
        findlowoff = find(lowertrans(t_offsets(i):tempend)<0,1,'first');
        if ~isempty(findlowoff)
            lowoffset(i) = findlowoff+t_offsets(i);
            tempstart = findlowoff+t_offsets(i);
        else
            lowoffset(i) = nan;
            tempstart = t_offsets(i);
        end

    end
    onsets = lowonset;
    offsets = lowoffset;

    %merge onsets and offsets without lower threshold crossings;
    onsets = onsets(~isnan(onsets));
    offsets = offsets(~isnan(offsets));
    assert(length(onsets)==length(offsets));
    assert(all((offsets-onsets)>0));
       
    %eliminate short intervals
    temp_int=(onsets(2:length(onsets))-offsets(1:length(offsets)-1))*1000/Fs;
    real_ints=temp_int>min_int;
    onsets=[onsets(1); nonzeros(onsets(2:length(onsets)).*real_ints)];
    offsets=[nonzeros(offsets(1:length(offsets)-1).*real_ints); offsets(length(offsets))];
    
    %eliminate short notes
    temp_dur=(offsets-onsets)*1000/Fs;
    real_durs=temp_dur>min_dur;
    onsets=[nonzeros((onsets).*real_durs)];
    offsets=[nonzeros((offsets).*real_durs)];
    
    %convert to ms: peculiarities here are to prevent rounding problem
    % if t_ons is simply replaced with onsets, everything gets rounded
    onsets = onsets/Fs; % all in seconds
    offsets = offsets/Fs; %all in seconds
end
return;
