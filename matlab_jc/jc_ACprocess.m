function del=jc_ACprocess(norm)
del=zeros(1,200);
for i=1:25
    covmat=xcov(norm(600:800,i));
    [peak,index]=max(covmat);
    for j=1:200
        del(j)=del(j)+(peak-covmat(index+j))/peak;
    end
end