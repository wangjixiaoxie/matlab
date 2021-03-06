%-1 get polygon bounds
disp(['Press return when doen with bounds\n']);

figure(CLUSTPLT);
if (get(DENPLOTBTN,'Value')==get(DENPLOTBTN,'Max'))
	set(DENPLOTBTN,'Value',get(DENPLOTBTN,'Min'));
	DensityPlot;
end

% get the polygon bounds
[x,y]=ginput;
% when return is pressed first point is added to close out the 
% bounds
x=[x;x(1)];y=[y;y(1)];
clust_hand=gca


if (length(x)<=3)
	%NO POLYGON
	return;
end

clr=get(COLORBOX,'String');
sym=get(SYMBOLBOX,'String');

hold on;
polyplt=plot(x,y,[clr,'--']);

pp=find(AMPPLOTS==gca);
chans = comb(pp,:);

insd=inpolygon(spkamp(:,chans(1)),spkamp(:,chans(2)),x,y);

%datnum is the fn number on the data plot
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
	
    %indices=find(spkind(pp,2)==177);
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

