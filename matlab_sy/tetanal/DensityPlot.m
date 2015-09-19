%function []=DensityPlot(a,b,c,CLUSTPLT,spkamp,DBINBOX);
%[]=DensityPlot(CLUSTPLT,spkamp,DBINBOX);
%
%
figure(CLUSTPLT);
if (get(DENPLOTBTN,'Value')==0)
	for ii = 1:size(comb,1)
		subplot(3,2,ii);hold off;
		plot(spkamp(:,comb(ii,1)),spkamp(:,comb(ii,2)),...
		                                    '.','MarkerSize',1);
		title([num2str(comb(ii,1)),' vs.' ,num2str(comb(ii,2))]);
		hold on;
	end
else
	if (isnumeric(str2num(get(DBINBOX,'String'))))
		dBIN = str2num(get(DBINBOX,'String'));
	else
		dBIN = 20;
		set(DBINBOX,'String',num2str(dBIN));
	end

	for ii=1:6
		DENSPLOTS(ii)=subplot(3,2,ii);hold off;
		inds = zeros([size(spkamp,1),2]);
		axbnds = zeros([2,1]);
		for jj = 1:2
			tmp=spkamp(:,comb(ii,jj));
			inds(:,jj)=floor((tmp-min(tmp))/dBIN)+1;
			axbnds(jj)=min(tmp);
		end
		denim = zeros([max(inds(:,2)),max(inds(:,1))]);
		for jj = 1:size(inds,1)
			denim(inds(jj,2),inds(jj,1))=...
			denim(inds(jj,2),inds(jj,1))+1;
		end
		imagesc([0:max(inds(:,2))-1]*dBIN+axbnds(1),...
		        [0:max(inds(:,1))-1]*dBIN+axbnds(2),denim);
		colormap('hot');
		m=colormap;m(1,:)=[1,1,1];
		colormap(m);set(gca,'YD','n');
		%v=axis;xlim([0,v(2)]);ylim([0,v(4)]);
		title([num2str(comb(ii,1)),' vs.' ,num2str(comb(ii,2))]);
	end
end
