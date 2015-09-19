function [Matrix]=jcnext(pitch,se)
wide=40;
for i=1:length(pitch)
    for j=1:length(se)
        Matrix(j,i)=median(pitch(i).pitches(se(j):se(j)+wide));
    end
end
