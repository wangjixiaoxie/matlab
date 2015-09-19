function jc_plotresidualstest731(norm,start)


sliced=norm(start,:);
xval(1)=1;
g(1)=std(sliced);



for j=1:15
    shift=40*j;     %5ms each time
    for i=1:length(sliced)
        dev(i)=norm(start+shift,i)-sliced(i);

    end
    xval(j+1)=shift;
    g(j+1)=std(dev);
end
hold on;
plot(xval,g,'MarkerSize',15,'Marker','.')
xlim([-40 600]);ylim([-0.05 0.03])
plot(norm(start:length(norm),1))