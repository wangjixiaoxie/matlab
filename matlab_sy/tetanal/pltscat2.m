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
    MKCLUSTBTN=uicontrol('Style','toggle','Position',[470,YP,70,20],...
				'String','Make Cluster','Value',0,...
				'Callback','clusterondensity');

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

if(clustnum>1)
%Separate spikes into two groups, those that are in the existing cluster,
%and those that are not in the existing cluster
%First we need to evaluate this with borrowe code from new cluster
    pp=find(AMPPLOTS==clust_hand);
    chans = comb(pp,:);
    insd=inpolygon(ampbuffer(:,chans(1)),ampbuffer(:,chans(2)),x,y);
%datnum is the fn number on the data plot
    ppin=find(insd==1);
    ppout=find(insd==0);
    inspikes=ampbuffer(ppin,:);
    outspikes=ampbuffer(ppout,:);

else
    
    inspikes=ampbuffer;
    outspikes=ampbuffer;
    

end


if (get(DENPLOTBTN,'Value')==0)
	for ii = 1:size(comb,1)
		
        scathandle(ii)=subplot(3,2,ii);
		%subplot(NCOL,ceil(size(comb,1)/NCOL),ii);
		plot(inspikes(:,comb(ii,1)),inspikes(:,comb(ii,2)),...
                                            'r.','MarkerSize',1);
                                        hold on;
		%Add a different color for the new trial
        plot(outspikes(:, comb(ii,1)),outspikes(:,comb(ii,2)),...
                                            'y.','MarkerSize',1);
                                            
        title([num2str(comb(ii,1)),' vs.' ,num2str(comb(ii,2))]);
    end
    %Add existing cluster if there is one.
    if(exist('clust_hand'))
        axes(clust_hand)
       polyplt=plot(x,y,['--']);
    end
else
	%DensityPlot(CLUSTPLT,spkamp,DBINBOX);
	DensityPlot
end

set(CLUSTPLT,'Toolbar','figure');
 
     