function [sumplot]=plotbarinitind3(sumplot, type,tag)

% if type='all'
%     
% elseif type='shi'
%         
% end
%first separate out upshifts and downshifts.
indup=find(sumplot(1).drxnvec>0)
inddn=find(sumplot(1).drxnvec==0)

if (exist('tag'))
    bstag=tag.bsnum
    colvls=tag.colvls;
    for cursp=1:length(sumplot)
        for vcnum=1:length(sumplot(cursp))
            sumplot(cursp).colvec{vcnum}='k';
        end
    end
    
    for ii=1:length(bstag)
        for cursp=1:length(sumplot)
            ind=find(sumplot(cursp).bsnum==bstag(ii));
            sumplot(cursp).colvec(ind)=colvls{ii}
        end
    end
end
    
    

vls{1}=sumplot(1).off(indup)
pctvls{1}=sumplot(1).pct(indup);
vls{2}=sumplot(1).off(inddn);
pctvls{2}=sumplot(1).pct(inddn);

for ii=1:2

    indlim=find(isnan(vls{ii})==0)
    mnoff{ii}=mean(vls{ii}(indlim))
    stdoff{ii}=std(vls{ii}(indlim))
    mnpct{ii}=mean(pctvls{ii}(indlim));
    stdpct{ii}=std(pctvls{ii}(indlim));
end

for ii=2:3
    if(ii==2)
        inind=ii;
        outind=3;
    else
        inind=4;
        outind=4;
    end
    vls{outind}=sumplot(inind).off;
    indlim=find(isnan(vls{outind})==0)
    mnoff{outind}=mean(vls{outind}(indlim))
    stdoff{outind}=std(vls{outind}(indlim))
end

if (type=='bar')
    range=[1:4]
    plotvls=vls;
    plotmean=mnoff;
    plotstd=stdoff;
else
    range=1:2
    plotvls=pctvls;
    plotmean=mnpct;
    plotstd=stdpct;
end
figure
for ii=range
   crvls=plotvls{ii};
   ln=length(crvls);
   xvals=ii*ones(ln,1);
   plot(ii,crvls,'ko');
   hold on;
   bar(ii,plotmean{ii},'FaceColor','none');
   plot([ii ii], [plotmean{ii}-plotstd{ii} plotmean{ii}+plotstd{ii}],'k');
end

if exist('tag')
    plot(plot_tag,plotvls{plot_tag}(tagindrvs),'ro')
end
 
        