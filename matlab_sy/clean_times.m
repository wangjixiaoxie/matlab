function [t_on_out,t_off_out] = clean_times(t_on,t_off,min_int,Fs);
% USAGE : [t_on_out,t_off_out] = clean_times(t_on,t_off,min_int,Fs);
% this is just a subroutine of the code in segment that
% cleans out the short intervals
temp_int =(t_on(2:length(t_on))-t_off(1:length(t_off)-1))*1000/Fs;
real_ints=(temp_int > min_int);
t_on_out =[t_on(1); nonzeros(t_on(2:length(t_on)).*real_ints)];
t_off_out=[nonzeros(t_off(1:length(t_off)-1).*real_ints);t_off(length(t_off))];
return;
