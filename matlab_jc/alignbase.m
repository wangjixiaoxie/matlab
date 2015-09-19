function alignbase(pd715)
figure; hold on;
for i=1:40
    k=pd715(i).pitches;
    for j=1:50
        s(j)=std(k(j*20:j*20+60));
    end
    plot(s)
end
