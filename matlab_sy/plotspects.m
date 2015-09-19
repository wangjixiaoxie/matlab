function plotspects(avn,t,f,lim,titles)

numspects=length(avn);
colormap('jet')


for i=1:numspects
    subplot(numspects,1,i)
    imagesc(t{i},f{i},log(avn{i}));syn;ylim([0,1e4]);
    grid on;
    box off;
    if(nargin>3)
        axis(lim);
        title(titles{i})
    end
    end