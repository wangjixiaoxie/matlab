function [t_on,t_off]=evsegment2(smooth,Fs,min_int,min_dur,max_dur,threshold);
% [t_on,t_off]=evsegment2(smooth,Fs,min_int,min_dur,max_dur,threshold);
% i was having trouble getting segemtn to find the quieter introductory
% sylables using a threshold that was high enough to avoid grouping the
% louder sylables together while low enough to catch these quieter notes
% here the threshold will be changed based on the local data statistics

disp(['BEWARE USING EVSEGMENT2!!!']);

%smooth the input file with a moving average filter of length min_dur
ndur = ceil(min_dur*Fs*1.0e-3);
sm2 = filter(ones([ndur,1]),1.0*ndur,smooth);

%get rid of convolution induced offset/filter delay effect
sm2=sm2(floor(ndur*0.5):end);

% get the first pass times of the threshold crossings using a
% low threhold value on the heavily filtered sm2 data set
n_tms = (sm2 > threshold);

%extract index values for note onsets and offsets
trans = diff(n_tms);
t_on  = find(trans > 0);
t_off = find(trans < 0);

% this should deal with a weird size mismatch
if (length(t_on) ~= length(t_off))
	disp('Number of note Onsets and Offsets do not match');
	if (length(t_on)<length(t_off)) 
		if (t_on(1)>t_off(1))
			t_off(1) = [];
		else
			disp(['Could not fix it!']);
			return;
		end
	else
		if (t_on(end)>t_off(end))
			t_on(end) = [];
		else
			disp(['Could not fix it!']);
			return;
		end
	end
	if (length(t_on) ~= length(t_off))
		disp(['Could not fix it!']);
		return;
	end
end
	
% the super smoothed sm2 should not have too many overly short intervals
% but get rid of them all the same
temp_int  = (t_on(2:length(t_on))-t_off(1:length(t_off)-1))*1000/Fs;
real_ints = (temp_int > min_int);
t_on  = [t_on(1); nonzeros(t_on(2:length(t_on)).*real_ints)];
t_off = [nonzeros(t_off(1:length(t_off)-1).*real_ints);t_off(length(t_off))];

size(t_on),size(t_off)

%start_point.mat

% some of these chunks are good others have a bunch of data in them
%reduce chunks > max_dur in length
temp_threshold_frac = 0.001;
durs = (t_off-t_on)*1000.0/Fs;
p = find(durs>=max_dur);
    
n_ton=[];n_toff=[];
for ii = 1:length(p)
    disp(num2str(ii));
    ind1 = t_on(p(ii));ind2 = t_off(p(ii));

    tdat = smooth([ind1:ind2]);
    tmp_thresh = (max(tdat)-min(tdat))*temp_threshold_frac + min(tdat);
    tmp_thresh = max([tmp_thresh,threshold]);

    tmp_ntimes = (tdat > tmp_thresh);
    tmp_trans = diff(tmp_ntimes);
	    
    new_off_times = find(tmp_trans < 0) + ind1 - 1;
    new_on_times  = find(tmp_trans > 0) + ind1 - 1;

    [new_on_times,new_off_times] = ...
         clean_times(new_on_times,new_off_times,min_int,Fs);

    n_ton  = [n_ton; new_on_times];
    n_toff = [n_toff;new_off_times];
    [n_ton,n_toff] = check_on_off_lengths(n_ton,n_toff);
end
    
t_on(p) = [];t_off(p) = [];
t_on = sort([t_on;n_ton]);t_off = sort([t_off;n_toff]);

size(t_on),size(t_off)
[t_on,t_off] = clean_times(t_on,t_off,min_int,Fs);
%eliminate short notes
temp_dur  = (t_off-t_on)*1000/Fs;
real_durs = (temp_dur > min_dur);
t_on  = [nonzeros((t_on).*real_durs)];
t_off = [nonzeros((t_off).*real_durs)];
    
%convert to ms: peculiarities here are to prevent rounding problem
% if t_ons is simply replaced with onsets, everything gets rounded
%t_on  = t_on*1000/Fs;
%t_off = t_off*1000/Fs;
return;
