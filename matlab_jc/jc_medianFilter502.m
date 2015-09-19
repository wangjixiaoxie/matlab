function [median_filtered]=jc_medianFilter502(input_matrix,halfwidth)
%use a halfwidth of 4 for optimal results.
width=halfwidth*2;
a=input_matrix;
c=a;



fudge_factor=80;
outlier=fudge_factor;

ending=length(a)-width;
for i=halfwidth+1:ending
    m1=median(a(i-halfwidth:i-1));
    m2=median(a(i+1:i+halfwidth));
    for j=0:halfwidth
        n=i+j;
        low1=m1-outlier;
        high1=m1+outlier;
        low2=m2-outlier;
        high2=m2+outlier;
        if (a(n)>high1 || a(n)<low1) && (a(n)>high2 || a(n)<low2)
            c(n)=(m1+m2)/2;
        end
    end
end
median_filtered=c;