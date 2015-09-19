function [median_filtered]=jc_medianFilter502(input_matrix,halfwidth)
%halfwidth of 3 gives a filter size of 6
width=halfwidth*2;
a=input_matrix;
c=a;
ending=length(a)-width;
for i=halfwidth+1:ending
    m=median(a(i-halfwidth:i+halfwidth));
    for j=0:width
        if a(i+j)>m+40 || a(i+j)<m-40
            c(i+j)=m;
        end
    end
end
median_filtered=c;