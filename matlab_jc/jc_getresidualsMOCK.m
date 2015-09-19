function [xcresids]=jc_getresidualsMOCK(processed,xmax,starting,ending)
%800,280,680 or 1000,400,800
note_length=6400;
sampling=32000;
ymin=0;
%Get the residuals
for i=1:length(processed(1).pitches)
    for j=1:length(processed)
        matr(i,j)=processed(j).pitches(i);
    end
    averaged(i)=mean(matr(i,:));
    for j=1:length(processed)
        normalized(i,j)=matr(i,j)/averaged(i);
        normalized(i,j)=normalized(i,j)-1;
    end
end

%Take the cross-correlation of the residuals for each trace.
for j=1:length(processed)
    xcr(:,j)=xcov(normalized(:,j));
end


for i=1:length(xcr)
    xcresids(i)=mean(xcr(i,:));
end