function [beginning,ending]=jc_segment(fvals)
fs=32000;
for i=1:length(fvals)
    ss=fvals(i).datt(1:(1024+round((fs/1000)*(fvals(i).offs(fvals(i).ind)-fvals(i).ons(fvals(i).ind)))));
    ending(i,:)=ss(length(ss)-2500:length(ss));
    beginning(i,:)=ss(1:2500);
end
    