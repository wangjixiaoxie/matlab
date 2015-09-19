function averaged=jc_plotresiduals725(pitchdata)



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

subplot(143);hold on;title('residuals');
for i=1:10
    alpha=normalized(:,i)+0.05*i;
    plot(alpha)
    plot(t+0.05*i,'LineStyle',':','Color',[1 0 0])
end
subplot(144);hold on;
for i=11:20
    alpha=normalized(:,i)+0.05*(i-10);
    plot(alpha)
    plot(t+0.05*(i-10),'LineStyle',':','Color',[1 0 0])
end