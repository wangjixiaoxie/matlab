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
	%set(COLORBOX,'String',COLORS(NCOLOR));
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
    for jj=1:totclusts
    
        pp(jj)=find(AMPPLOTS==clust_hand{jj});
        chans{jj} = comb(pp(jj),:);
        chans2=chans{jj}
        
        insd{jj}=inpolygon(ampbuffer(:,chans2(1)),ampbuffer(:,chans2(2)),x{jj},y{jj});
%datnum is the fn number on the data plot
        ppin{jj}=find(insd{jj}==1);
        ppout{jj}=find(insd{jj}==0);
        inspikes{jj}=ampbuffer(ppin{jj},:);
        outspikes=ampbuffer;
    end
else
    %This is just to make inspikes a cell array
    inspikes{1}=ampbuffer;
    outspikes=ampbuffer;
    

end


%if (get(DENPLOTBTN,'Value')==0)
	for ii = 1:size(comb,1)
		scathandle(ii)=subplot(3,2,ii);
        hold on;
        plot(outspikes(:, comb(ii,1)),outspikes(:,comb(ii,2)),...
                                            'y.','MarkerSize',1);
        for jj=1:totclusts
            
            %subplot(NCOL,ceil(size(comb,1)/NCOL),ii);
           
            
            if(clustnum>1)
                plot(inspikes{jj}(:,comb(ii,1)),inspikes{jj}(:,comb(ii,2)),...
                                            '.','MarkerSize',1,'Color',colorlist{jj});
                
            end                           
                                            
            
            %Add a different color for the new trial
            title([num2str(comb(ii,1)),' vs.' ,num2str(comb(ii,2))]);
        
    %Add existing cluster if there is one.
            if(skiptrial~=0)
                axes(clust_hand{jj})
                polyplt=plot(x,y,['--']);
            end
        end
    end



            %else



	%DensityPlot(CLUSTPLT,spkamp,DBINBOX);
%	DensityPlot
%end

set(CLUSTPLT,'Toolbar','figure');
              
     