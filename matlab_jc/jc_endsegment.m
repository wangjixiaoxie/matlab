function [ss]=jc_endsegment(fvals,sig,first,last)
fs=32000;
for i=1:length(fvals)
    ss=fvals(i).datt(1:(1024+round((fs/1000)*(fvals(i).offs(fvals(i).ind)-fvals(i).ons(fvals(i).ind)))));
    ending(i,:)=ss(length(ss)-2500:length(ss));
    beginning(i,:)=ss(1:2500);
end
%pitch=jc_PitchData729_1024(ending,1024,1016,sig,first,last);
%for i=1:length(pitch)
%    aa(i)=median(pitch(i).pitches(140:180));
%end
    
    