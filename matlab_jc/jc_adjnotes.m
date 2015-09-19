function [m,n,nn]=jc_adjnotes(norm)
j=1;
for i=1:size(norm,2)
    m(i)=median(norm(1800:2000,i));

    n(i)=median(norm(2600:2800,i));
    if m(i)>-0.008
        nn(j)=n(i);
        j=j+1;
    end
    mn(i)=n(i)-m(i); %distance between them
end
%guess: the std of mn is less than the std of n.