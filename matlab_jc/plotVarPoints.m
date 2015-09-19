function plotVarPoints(Experiment7,VarPre,VarPost,VarAPV,VarPost1,VarPost2,VarPreCTL,VarPostCTL,VarAPVCTL,VarPost1CTL,VarPost2CTL,ind,ind2)

        figure;hold on;
     VarPost2CTL(7)=83.1951/mean(mean(Experiment7.pitchACpostCTL(170:340,:))); % after removing outliers   
        subplot(121);hold on;plot(VarPre(ind),VarAPV(ind),'r.','Markersize',15);plot(VarPre(ind),VarPost2(ind),'b.','Markersize',15);plot([0 0.04],[0 0.04],'k');%plot(VarPre(ind),VarPost1(ind),'g.');
 mepre=mean(VarPre(ind));
 sepre=std(VarPre(ind))/sqrt(length(VarPre(ind)));
  meapv=mean(VarAPV(ind));
 seapv=std(VarAPV(ind))/sqrt(length(VarAPV(ind)));
 mepost=mean(VarPost2(ind));
 sepost=std(VarPost2(ind))/sqrt(length(VarPost2(ind)));
        plot([mepre mepre],[meapv-seapv meapv+seapv],'r')
        plot([mepre-sepre mepre+sepre],[meapv meapv],'r')
        plot([mepre mepre],[mepost-sepost mepost+sepost],'b')
        plot([mepre-sepre mepre+sepre],[mepost mepost],'b')

        
        
        subplot(122);hold on;plot(VarPreCTL(ind),VarAPVCTL(ind),'r.','Markersize',15);plot(VarPreCTL(ind),VarPost2CTL(ind),'b.','Markersize',15);plot([0 0.05],[0 0.05],'k')   
 mepreCTL=mean(VarPreCTL(find(VarPreCTL>0)));
 sepreCTL=std(VarPreCTL(find(VarPreCTL>0)))/sqrt(length(VarPreCTL(find(VarPreCTL>0))));
  meapvCTL=mean(VarAPVCTL(find(VarAPVCTL>0)));
 seapvCTL=std(VarAPVCTL(find(VarAPVCTL>0)))/sqrt(length(VarAPVCTL(find(VarAPVCTL>0))));
 mepostCTL=mean(VarPost2CTL(find(VarPost2CTL>0)));
 sepostCTL=std(VarPost2CTL(find(VarPost2CTL>0)))/sqrt(length(VarPost2CTL(find(VarPost2CTL>0))));
        plot([mepreCTL mepreCTL],[meapvCTL-seapvCTL meapvCTL+seapvCTL],'r')
        plot([mepreCTL-sepreCTL mepreCTL+sepreCTL],[meapvCTL meapvCTL],'r')
        plot([mepreCTL mepreCTL],[mepostCTL-sepostCTL mepostCTL+sepostCTL],'b')
        plot([mepreCTL-sepreCTL mepreCTL+sepreCTL],[mepostCTL mepostCTL],'b')
