%stability_analysis
%takes input from ~/ml/seqbirdmfiles/synbirdstruct
function [outmat,ss]=stabil_plot(sumout,ss,PLOTFLAG)



%first plot example cbin
if(PLOTFLAG)
    subplot(3,3,1:3)
    ax=gca();
    exsong.ax=ax;
    exsong.path='/oriole7/dir2/pu19bk81/screen'
    exsong.fn='pu19bk81_280311_1254.-29188.cbin';
    exsong.bnds=[3.75 5.65]
    plotcbin(exsong);

    %this plots example transition probabilities for 6 days for this example
    subplot(3,3,4:5);
    exind=6;
    crdt=sumout(exind);
    ps.days=crdt.days;
    ndys=length(crdt.days)
    ps.prob(:,1)=[100*ones(1,ndys)*(crdt.sqvls)]
    ps.prob(:,2)=[100*ones(1,ndys)*(1-crdt.sqvls)]
    ps.colors={'k' 'r'}
    ps.err=[];
    axis([0 7 0 100])
    set(gca,'ytick',[0:20:100])
end
dyinds=[1:4]
for ii=1:length(sumout)
   
    crdt=sumout(ii);

        inds=find(ismember(crdt.days, dyinds));
        mnout=mean(crdt.sqvls(inds));
       absdiffvls(ii,:)=abs(((crdt.sqvls(inds)/mnout)-1)*100);
    if(PLOTFLAG)
       subplot(3,3,9)
       plot(crdt.sqvls(dyinds(1))*100,crdt.sqvls(dyinds(4))*100,'ko');
       hold on;
    end
    outmat.ssind(ii)=sumout(ii).ssind;   
    ss(sumout(ii).ssind).basprob=mnout;
end
if(PLOTFLAG)
    plot([0 100],[0 100],'k--')
    set(gca,'ytick',[0:20:100]);
    set(gca,'xtick',[0:20:100]);
    axis square;
    box off;
    subplot(3,3,7:8);

        mnvls_sum=mean(absdiffvls,1);
        ste_sum=std(absdiffvls,0,1)./sqrt(length(sumout));
        bar([dyinds],mnvls_sum,0.6,'Facecolor','none');
        hold on;
        plot([dyinds;dyinds],[(mnvls_sum+ste_sum);(mnvls_sum-ste_sum)],'k');
    axis([0 5 0 20])
    set(gca,'ytick',[0:5:20])
    hold on;
    plot([1 4],[mean(mnvls_sum) mean(mnvls_sum)],'k--');
    box off;
    axis square;
end
outmat.vls=absdiffvls;


%         plot([1 2 3],1-mnout,'r');
%         hold on;
%         plot([1 2 3],1-mnout,'ro');
%         
%         subplot(3,3,7);
%         plot([1 2 3],mean(mnout)-mnout,'k');
%         hold on;
%         plot([1 2 3],mean(mnout)-mnout,'ko');

%         plot([1 2 3],1-(mean(mnout)-mnout),'r');
%         hold on;
%         plot([1 2 3],1-(mean(mnout)-mnout),'ro');

   
