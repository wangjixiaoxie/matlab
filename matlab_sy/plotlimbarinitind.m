function [sumplot]=plotlimbarinitind(sumplot, type)

% if type='all'
%     
% elseif type='shi'
%         
% end

i

%first separate out upshifts and downshifts.
indup=find(sumplot(1).ac>0)
inddn=find(sumplot(1).ac<0)

indup=intersect(indup,sumplot(1).indthresh);
inddn=intersect(inddn,sumplot(1).indthresh);

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
    vls{ii+1}=sumplot(ii).off;
    indlim=find(isnan(vls{ii+1})==0)
    mnoff{ii+1}=mean(vls{ii+1}(indlim))
    stdoff{ii+1}=std(vls{ii+1}(indlim))
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
 
        