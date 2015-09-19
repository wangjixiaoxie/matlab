function [pitch,shifted]=pitchcontour(fvals,minwin,maxwin)
% pitch=pitchcontour(fvals,2000,3000);

for i=1:length(fvals)
    shifted(i,:)=fvals(i).datt;
end
pitch=jc_pitchmat1024(shifted,1024,1020,1,minwin,maxwin,[1],'obs0',1);
