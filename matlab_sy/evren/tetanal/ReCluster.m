figure(CLUSTPLT);
if (get(DENPLOTBTN,'Value')==get(DENPLOTBTN,'Max'))
	set(DENPLOTBTN,'Value',get(DENPLOTBTN,'Min'));
	DensityPlot;
end

CLUSTNUM = str2num(get(CURRCLUSTERBOX,'String'));
x=clustdat(CLUSTNUM).polygon(:,1);
y=clustdat(CLUSTNUM).polygon(:,2);
chans=clustdat(CLUSTNUM).polygonchans;

clr=clustdat(CLUSTNUM).clr;
sym=clustdat(CLUSTNUM).sym;

if (length(x)<=3)
	%NO POLYGON
	return;
end

pp=find((comb(:,1)==chans(1))&(comb(:,2)==chans(2)));
chans = comb(pp,:);

subplot(3,2,pp);
hold on;
polyplt=plot(x,y,[clr,'--']);

insd=inpolygon(spkamp(:,chans(1)),spkamp(:,chans(2)),x,y);

pp=find(insd==1);
for ii = 1:size(comb,1)
	%subplot(NCOL,ceil(size(comb,1)/NCOL),ii);
	subplot(3,2,ii);
	hold on;
	plot(spkamp(pp,comb(ii,1)),spkamp(pp,comb(ii,2)),...
	            [clr,'.'],'MarkerSize',1);
end

figure(DATPLT);
for ii = 1:Nusechans
	subplot(Nusechans+1,1,ii+1);
	hold on;grid on;
	%plot(spkt(pp),dat(spki,usechans(ii-1)),[clr,'o']);
	plot(spkind(pp)/fs,spkamp(pp,ii),[clr,sym]);
end

set(CURRCLUSTERBOX,'String',num2str(CLUSTNUM+1));

NCOLOR=CLUSTNUM+1;
if (CLUSTNUM<length(COLORS))
	set(COLORBOX,'String',COLORS(NCOLOR));
end