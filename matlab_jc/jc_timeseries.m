function [x,y]=jc_timeseries(medians)
for i=5:length(medians)
    x(i-4)=medians(i-4);
    y(i-4)=medians(i);
end
figure; plot(x,y,'*')
