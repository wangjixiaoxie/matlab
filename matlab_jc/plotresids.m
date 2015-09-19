function plotresids(norm)

subplot(131); xlabel('Time (ms)'); xlim([0 800]);hold on;
for i=1:10
    norm(:,i)=norm(:,i)+0.03*i;
    plot(norm(:,i))
    
end

subplot(132); xlabel('Time (ms)'); xlim([0 800]);hold on;
for i=11:20
    norm(:,i)=norm(:,i)+0.03*(i-10);
    plot(norm(:,i))
end

subplot(133); xlabel('Time (ms)'); xlim([0 800]);hold on;
for i=21:30
        norm(:,i)=norm(:,i)+0.03*(i-20);
        plot(norm(:,i))
end