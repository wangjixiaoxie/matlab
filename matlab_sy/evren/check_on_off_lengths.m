function [t_on,t_off] = check_on_off_lengths(t_on_tmp,t_off_tmp);
% [t_on,t_off] = check_on_off_lengths(t_on_tmp,t_off);
% correct if the length of t_on_tmp and t_off_tmp are not the same
% 
if (length(t_on_tmp) ~= length(t_off_tmp))
	disp('Number of note Onsets and Offsets do not match');
	if (t_on_tmp(1)>t_off_tmp(1))
		t_off_tmp(1) = [];
	end

	if (t_on_tmp(end)>t_off_tmp(end))
		t_on_tmp(end) = [];
	end

	if (length(t_on_tmp) ~= length(t_off_tmp))
		disp(['Could not fix it!']);
		t_on = t_on_tmp;t_off = t_off_tmp;
		return;
	end
end
t_on = t_on_tmp;t_off = t_off_tmp;
return;
