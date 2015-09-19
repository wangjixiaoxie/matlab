function ave=jcwonders717b(fv1)
figure;
hold on;
b=zeros(6500,1);
j=0;
for i=1:length(fv1)
    a=fv1(i).datt;
    a=SmoothData(a,32000,1);
    a=log(a);
    plot(a)
    golf(i).data=a;
    if(a(4000)<18)
        b=b+a;
        j=j+1;
    end
end
ave=b/j;
plot(ave)
