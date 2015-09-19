% figure
% ax(1)=subplot(212)
% for ii=1:length(vlsout)
% plot(mean(vlsout{ii}(:,1)),mean(vlsout{ii}(:,2)),'ro','Markersize',6)
% hold on;
% end
% 
figure
ax(2)=subplot(211)
for ii=1:length(vlsout)
plot(mean(tms{ii}),mean(pout{ii}),'ro','MarkerSize',6)
hold on;
end
linkaxes(ax,'x')

figure
for ii=1:length(vlsout)
xvls=[tms{ii};tms{ii}(end:-1:1)]
yvls=[p05out{ii},p95out{ii}(end:-1:1)]
fill(xvls,yvls,'r')
hold on;
end

