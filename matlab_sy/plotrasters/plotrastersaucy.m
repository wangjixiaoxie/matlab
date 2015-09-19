figure

sn=1
nfiles=6

%choose files to plot
outind=rand(nfiles,1);
pltind=ceil(outind*stimf(sn).cnt)

for ii=1:length(pltind)
    ax6(ii)=subplot(nfiles+1,1,ii)
    dat=ReadCbinFile(fnames(pltind(ii)).fn);
    plot(0:1/fs:(length(dat(:,3))-1)/fs,dat(:,3))
    hold on;
    plotspkarray(stimf(sn).spkarray{pltind(ii)}*32000, -10000, 1000);
        
    title(fnames(pltind(ii)).fn)
end    
linkaxes(ax6);

ax7=subplot(7,1,7)
plot(0:1/fs:(length(dat(:,3))-1)/fs,dat(:,1))


