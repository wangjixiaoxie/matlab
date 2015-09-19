function jc_timeseries(medians)
for i=2:length(medians)
    x=medians(i-1);
    y=medians(i);
end
figure; plot(x,y)
    