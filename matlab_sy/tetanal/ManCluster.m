%-1 get polygon bounds

figure(CLUSTPLT);
if (get(DENPLOTBTN,'Value')==get(DENPLOTBTN,'Max'))
	set(DENPLOTBTN,'Value',get(DENPLOTBTN,'Min'));
	DensityPlot;
end

clr=get(COLORBOX,'String');
sym=get(SYMBOLBOX,'String');

pp=find(AMPPLOTS==gca);
chans = comb(pp,:);

insd=inpolygon(spkamp(:,chans(1)),spkamp(:,chans(2)),x,y);

pp=find(insd==1);
for ii = 1:size(comb,1)
	subplot(3,2,ii);
	hold on;
	plot(spkamp(pp,comb(ii,1)),spkamp(pp,comb(ii,2)),...
	            [clr,'.'],'MarkerSize',1);
end

figure(DATPLT);
for ii = 1:Nusechans
	subplot(Nusechans+1,1,ii+1);
	hold on;grid on;
	plot(spkind(pp)/fs,spkamp(pp,ii),[clr,sym]);
end


%SAVE ALL THE INFO
NCLUST=NCLUST+1;
clustdat(NCLUST).sym = sym;
clustdat(NCLUST).clr = clr;
clustdat(NCLUST).polygon = [x,y];
clustdat(NCLUST).polygonchans = chans;
clustdat(NCLUST).spikes = pp;

set(CURRCLUSTERBOX,'String',num2str(NCLUST));

if (NCOLOR<length(COLORS))
	NCOLOR=NCOLOR+1;
	set(COLORBOX,'String',COLORS(NCOLOR));
end

