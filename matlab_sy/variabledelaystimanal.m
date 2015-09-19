fvcatchcontour=[]
fvstimtime=[];
fvdel=[];

for ii=1:length(fv)
    if((fv(ii).STIMTRIG==1))
        fvstimtime=[fvstimtime fv(ii).STIMTIME]
     if(fv(ii).STIMCATCH==1)
        fvcatchcontour=[fvcatchcontour;fv(ii).pitch_data];
     else
         fvdel=[fvdel fv(ii).STIMTIME]
     end
             
         
    end
end

make residual
for ii=1:length(fv)
    fv(ii).resid=mn-fv(ii).pitch_data;
end
edges=[46:5:76]
ln=length(edges)
figure
for ii=2:length(edges)
    ind=find(fvstimtime<edges(ii)&fvstimtime>edges(ii-1))
    pdcomb{ii}=[];
    for jj=1:length(ind)
        crind=ind(jj);
        pdcomb{ii}=[pdcomb{ii};fv(ii).pitch_data];
    end
    mnpd{ii}=mean(pdcomb{ii},1);
    ax(ii-1)=subplot(ln-1,1,ii-1)
    if(~isempty(mnpd{ii}))
        plot(pitchtms,mn-mnpd{ii})
    end
end
    
end

for 
    
