function  shadeplot(xdat, datmtx,colidx)
load lenacols2
colsi = cols;
load ctxtcols
colsi(5:8,:) = cols;
cols = colsi;

% load errorcols
% datacol = [.2 .2 1];
shadecol = [.6 .8 1; .6 .8 1; 1 .8 .6; 1 .8 .6];


ndata = sum(~isnan(datmtx),1);

xdat(ndata==0) = [];
datmtx(:,ndata==0) = [];
ndata(ndata==0) = [];

dat = nanmean(datmtx,1);
semi = nanstd(datmtx,1)./sqrt(ndata);

% dat = dat(~isnan(semi));
% semi = semi(~isnan(semi));
% xdat = xdat(~isnan(semi));

interval_lo = dat-semi;
interval_hi = dat+semi;

% h = patch([xdat,xdat(end:-1:1)],[interval_lo interval_hi(end:-1:1)],shadecol(colidx,:),'EdgeColor','none','FaceAlpha',.5)
hold on
for i = 1:size(dat,1)
%     plot(xdat,(dat(i,:)),'Color',cols(colidx,:),'LineWidth',1)
%     plot(xdat,(dat(i,:)),'Color',datacol,'LineWidth',3)

end
plot(xdat,dat,'-','color',cols(colidx,:),'linewidth',2)
h = errorbar(xdat,dat,semi,'color',cols(colidx,:),'linestyle','none');
set(get(get(h,'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
