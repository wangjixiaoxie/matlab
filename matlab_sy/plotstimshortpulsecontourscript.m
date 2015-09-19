figure
plotinds=[3 4  5 6]
stimt=[.030 .040 .050 .060]
ln=length(plotinds);
for ii=1:ln
    crind=plotinds(ii);
    cr_time=stimt(ii);
    ax(ii)=subplot(ln,1,ii)
    mnct=mean(outstr.combctcon{crind},2)
    mnfb=mean(outstr.combfbcon{crind},2)
    sterct=std(outstr.combctcon{crind},0,2)./sqrt(length(outstr.combctcon{crind}(1,:)));
    sterfb=std(outstr.combfbcon{crind},0,2)./sqrt(length(outstr.combfbcon{crind}(1,:)));
    plot(pitchtms,mnct,'k','Linewidth',3)
    hold on;
    plot(cr_time+stimtms(45:390)-stimtms(45),stimblockscale(45:390),'k','Linewidth',2)
    plot(pitchtms,mnct+sterct,'k','Linewidth',1)
    plot(pitchtms,mnct-sterct,'k','Linewidth',1)
    plot(pitchtms,mnfb,'r','Linewidth',3)
    hold on;
    plot(pitchtms,mnfb+sterfb,'r','Linewidth',1)
    plot(pitchtms,mnfb-sterfb,'r','Linewidth',1)
    plot(pitchtms,mnctbas,'k--','Linewidth',2);
end
 linkaxes(ax)   