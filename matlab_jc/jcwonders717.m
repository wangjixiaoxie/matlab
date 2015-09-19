function ave=jcwonders717(fv1)
figure;
hold on;
b=zeros(6500,1);
j=0;
for i=1:20
    a=fv1(i).datt;
    a=SmoothData(a,44100,1);
    a=log(a);
    golf(i).data=a;
    if(a(4000)<18)
        b=b+a;
        j=j+1;
    end
end
ave=b/j;

for i=1:20
    h=xcov(ave(1:3000),golf(i).data(1:3000));
    [peak,index]=max(h);
    shift(i)=(index-6500); %6400 is note length
    k=1;
    for j=1:5000
        align_shift=shift(i);
        shift_index(i)=round(j-align_shift);
        if shift_index(i)>0
            shifted(i,k)=golf(i).data(shift_index(i));
        else
            shifted(i,k)=100;
        end
        k=k+1;
    end
    plot(shifted(i,:));
end
plot(ave);
