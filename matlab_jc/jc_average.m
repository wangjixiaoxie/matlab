function [averaged,up,dn]=jc_average(pitchdata)



%Get the residuals
for i=1:length(pitchdata(1).pitches)
    t(i)=0;
    for j=1:length(pitchdata)
        matr(i,j)=pitchdata(j).pitches(i);
    end
    standev(i)=std(matr(i,:));
    averaged(i)=mean(matr(i,:));
    up(i)=averaged(i)+standev(i)/sqrt(length(pitchdata));
    dn(i)=averaged(i)-standev(i)/sqrt(length(pitchdata));
end

