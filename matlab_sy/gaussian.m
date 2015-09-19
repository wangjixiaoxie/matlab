%mock gaussian plot
figure
subplot(311)
x=-20:.0001:20
sig1=4;
sig2=2.5;
ave1=0;
y=1/(sqrt(2*pi*sig1^2))*exp(-(x-ave1).^2./(2*sig1.^2));
plot(x,y,'k','Linewidth',2);
hold off

hold on;
y=1/(sqrt(2*pi*sig2^2))*exp(-(x-ave1).^2./(2*sig2.^2));
plot(x,y,'r--','Linewidth',2);
box off;
plotvertline([0 0],[0 0.2],'--');

subplot(312)
x=-20:.0001:20
sig1=4;
sig2=2.5;
ave1=12;
y=1/(sqrt(2*pi*sig1^2))*exp(-(x-ave1).^2./(2*sig1.^2));
plot(x,y,'k','Linewidth',2);
hold off

hold on;
y=1/(sqrt(2*pi*sig2^2))*exp(-(x-ave1).^2./(2*sig2.^2));
plot(x,y,'r--','Linewidth',2);
box off;
plotvertline([0 0],[0 0.2],'--');


subplot(313)
x=-20:.0001:20
sig1=4;
sig2=2.5;
ave1=12;
ave2=6;
y=1/(sqrt(2*pi*sig1^2))*exp(-(x-ave1).^2./(2*sig1.^2));
plot(x,y,'k','Linewidth',2);
hold off

hold on;
y=1/(sqrt(2*pi*sig2^2))*exp(-(x-ave2).^2./(2*sig2.^2));
plot(x,y,'r--','Linewidth',2);
box off;
plotvertline([0 0],[0 0.2],'--');