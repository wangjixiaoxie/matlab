function normalized=jc_plotresiduals725a(pitchdata)



%Get the residuals
for i=1:length(pitchdata(1).pitches)
    t(i)=0;
    for j=1:length(pitchdata)
        matr(i,j)=pitchdata(j).pitches(i);
    end
    averaged(i)=mean(matr(i,:));
end

for i=1:length(pitchdata)
    for j=1:length(pitchdata(1).pitches)
        normalized(j,i)=matr(j,i)/averaged(j);
        normalized(j,i)=normalized(j,i)-1;
    end
end


figure;hold on;title('residuals');
for i=1:5
    alpha=normalized(:,i)+0.1*i;
    plot(alpha)
    plot(t+0.1*i,'LineStyle',':','Color',[1 0 0])
end
figure;hold on;
for i=6:10
    alpha=normalized(:,i)+0.1*(i-5);
    plot(alpha)
    plot(t+0.1*(i-5),'LineStyle',':','Color',[1 0 0])
end
figure;hold on;
for i=11:13
    alpha=normalized(:,i)+0.1*(i-10);
    plot(alpha)
    plot(t+0.1*(i-10),'LineStyle',':','Color',[1 0 0])
end