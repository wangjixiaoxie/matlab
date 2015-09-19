%rewrite 4.17.09 with option to combine upshifts and downshifts.
%1. make option to only select first inactivation run.
%2. select controls to match (obligatory)

%input of sumplot, 
%first field of struct is pitchshifted offset 
%second field is control
%third is baseline.
%fourth is summary values by birds.
function [vls]=plotbarinitind6(sumplot, type,sel_init,tag)

% if type='all'
%     
% elseif type='shi'
%         
% end
%first separate out upshifts and downshifts.
indupin=find(sumplot(1).drxnvec>0)
inddnin=find(sumplot(1).drxnvec==0)
if (exist('sel_init'))
    if(sel_init)
        indup=redind(sumplot(1),indupin);
        inddn=redind(sumplot(1),inddnin);
    else
        indup=indupin;
        inddn=inddnin
    end
else
    indup=indupin;
    inddn=inddnin;
end


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
 

%At this point, the first argument of struct is is all the target notes
%second is all the control notes.
%third is all the baseline notes.
%first output of vals is upshift target
%second output is downshifted target
%3rd is control upshift
%4th is control downshift
%5th is summary.
 
vls(1).off=sumplot(1).off(indup)
vls(1).pct=sumplot(1).pct(indup);
vls(1).bsnum=sumplot(1).bsnum(indup)
vls(2).off=sumplot(1).off(inddn);
vls(2).pct=sumplot(1).pct(inddn);
vls(2).bsnum=sumplot(1).bsnum(inddn);
vls(3).off=sumplot(2).off(indup)
% pctvls{3}=sumplot(2).pct(indup);
vls(4).off=sumplot(2).off(inddn);
% pctvls{4}=sumplot(2).pct(inddn);
for ii=1:4

    indlim=find(isnan(vls(ii).off)==0)
    mnoff{ii}=mean(vls(ii).off(indlim))
    stdoff{ii}=std(vls(ii).off(indlim))
    if (ii<3)
    mnpct{ii}=mean(vls(ii).pct(indlim));
    stdpct{ii}=std(vls(ii).pct(indlim));
    end
end

    for ii=3:3
    %shift summary baseline to 5th out.
        inind=4;
        outind=5;
        vls(outind).off=sumplot(inind).off;
        indlim=find(isnan(vls(outind).off)==0)
        mnoff{outind}=mean(vls(outind).off(indlim))
        stdoff{outind}=std(vls(outind).off(indlim))
    end

   


if (type=='bar')
    range=[1:5]
    plotvls=vls;
    plotmean=mnoff;
    plotstd=stdoff;
else
    range=1:2
    plotvls=vls;
    plotmean=mnpct;
    plotstd=stdpct;
end
figure
for ii=range
   if(type=='bar')
       crvls=plotvls(ii).off;
   else
       crvls=plotvls(ii).pct
   end
       
       ln=length(crvls);
   xvals=ii*ones(ln,1);
if type=='bar'
    xvals=[1 4 2 5 7]
else
    xvals=[1 2]
end
   plot(xvals(ii),crvls,'ko');
   hold on;
   bar(xvals(ii),plotmean{ii},'FaceColor','none');
   plot([xvals(ii) xvals(ii)], [plotmean{ii}-plotstd{ii} plotmean{ii}+plotstd{ii}],'k');
end

if exist('tag')
    plot(plot_tag,plotvls{plot_tag}(tagindrvs),'ro')
end
 
%purpose of this function is to pull out multiple inactivations in same
%shift, stay with the first one.
function [indinit]=redind(sumplot, ind)
    indinit=ind(1);
    for ii=2:length(ind)
        indvl=ind(ii)
        if ((sumplot.bsnum(indvl)==sumplot.bsnum(indvl-1))&(sumplot.shftnum(indvl)==sumplot.shftnum(indvl-1))&(sumplot.drxnvec(indvl)==sumplot.drxnvec(indvl-1)))
        else
            indinit=[indinit indvl]
        end
    end
