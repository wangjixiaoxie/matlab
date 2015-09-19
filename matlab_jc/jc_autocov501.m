function jc_autocov501(input_array,first_note,last_note,start_bin,stop_bin)
y_min=first_note*1.5-1;
for i=first_note:last_note
    a=input_array(i).pitches(start_bin:stop_bin);
    b=xcov(a)/max(xcov(a));  %Normalize the autocovariance matrix
    shifter=i*1.5;
    y_max=shifter+1;
    c=b+shifter;  %Shift the graph up
    plot(c); ylim([y_min y_max])
end
