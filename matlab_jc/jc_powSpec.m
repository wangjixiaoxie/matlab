function del=jc_powSpec(norm)
for i=1:25
    covmat=xcorr(norm(150:500,i));
    h(i,:)=jc_fft(covmat);
end
del=zeros(1,size(h,2));
for i=1:25
    for j=1:length(del)
        del(j)=del(j)+h(i,j);
    end
end
    