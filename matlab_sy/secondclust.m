clear clust_hand
display_skip=10;
spksinclust=[];
ampbuffer=[];
indbuffer=[];
offset=0;
clustnum=1;

for ifn=0:5:length(fnm)
    
        zerovec=zeros(length(spkind),1);
        zerovec(spksin{1},1)=1;
    
        indices=find(spkind(:,2)>ifn&spkind(:,2)<(ifn+display_skip));
    
    
        ampbuffer=[spkamp(indices,:)];
        zerovec=zerovec(indices,:);
    
    

    
        pltdat;
        pltscat3;
        
        while(1)
                clustnum=clustnum+1;
                R=input('Keep Cluster? 1 for Yes')
                %strname is the correct name for clustered spikes%
                if (R==1)
                    spksinclust=[spksinclust;ppin+offset];
                    %spksoutclust=[spksoutclust;ppout+offset];
                     totalpolygons=totalpolygons+1;
                     clusterinfo{totalpolygons}.dims=[x,y];
                     clusterinfo{totalpolygons}.chans=[chans];
                break;
                else
                    figure(CLUSTPLT);
                    [x,y]=ginput();
                    clust_hand=gca;
                    pltdat;
                    pltscat3;
                    
                    figure(DATPLT);
                        %for ii = 1:Nusechans
                         %   indices=find(inspikes(:,5)==ifn)
                          %  indices2=find(indbuffer(ppin,2)==ifn)
                           % subplot(Nusechans+1,1,ii+1);
                            %hold on;grid on;
                        
                        %end
             
                end
     
        end
    offset=offset+length(ampbuffer(:,1));
    ampbuffer=[];
    indbuffer=[];
    clustnum=clustnum+1;  
    end

