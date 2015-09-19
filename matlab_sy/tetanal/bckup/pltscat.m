comb = [1,2;1,3;1,4;2,3;2,4;3,4];

% cluster colors
NCOLOR=1;
COLORS=['rgmkcy'];

% save the cluster data
NCLUST=0;
clustdat=[];

if (~exist('CLUSTPLT'))
	figure;
	CLUSTPLT=gcf;
	set(CLUSTPLT,'Tag','ClusterPlot','Position',[676,464,595,485],...
	             'DeleteFcn','CLUSTPLOTDELETE');
	YP=5;
	SYMBOLBOX=uicontrol('Style','edit','Position',[70,YP,20,20],...
	            'String','o');
	COLORBOX=uicontrol('Style','edit','Position',[100,YP,20,20],...
	            'String','r');
	MAKENEWCLUST=uicontrol('Style','pushbutton',...
	                       'Position',[130,YP,70,20],...
			       'String','Cluster','Callback','NewCluster');
	DBINBOX=uicontrol('Style','edit','Position',[290,YP,30,20],...
				'String','100');
	DENPLOTBTN=uicontrol('Style','toggle','Position',[210,YP,70,20],...
				'String','Density','Value',0,...
				'Callback','DensityPlot');
	AVSPKBTN=uicontrol('Style','pushbutton','Position',[330,YP,70,20],...
			'String','AvSpk','CallBack','AvSpk');
	CURRCLUSTERBOX=uicontrol('Style','edit','Position',[410,YP,30,20],...
			'String','0');
else
	figure(CLUSTPLT);
	for ii = 1:6
		subplot(3,2,ii);cla;
	end
	set(COLORBOX,'String',COLORS(NCOLOR));
end

if (exist('AVSPKPLT'))
	 figure(AVSPKPLT);
	 for ii=1:4
		 subplot(2,2,ii);
		 cla;
	 end
	 figure(CLUSTPLT);
end

AMPPLOTS=zeros([6,1]);
for ii = 1:6
	AMPPLOTS(ii)=subplot(3,2,ii);
end

if (get(DENPLOTBTN,'Value')==0)
	for ii = 1:size(comb,1)
		subplot(3,2,ii);
		%subplot(NCOL,ceil(size(comb,1)/NCOL),ii);
		plot(spkamp(:,comb(ii,1)),spkamp(:,comb(ii,2)),...
		                                    '.','MarkerSize',1);
		title([num2str(comb(ii,1)),' vs.' ,num2str(comb(ii,2))]);
        end
else
	%DensityPlot(CLUSTPLT,spkamp,DBINBOX);
	DensityPlot
end

set(CLUSTPLT,'Toolbar','figure');
