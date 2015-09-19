function plotVarLines(Experiment7,VarPre,VarPost,VarAPV,VarPost1,VarPost2,VarPreCTL,VarPostCTL,VarAPVCTL,VarPost1CTL,VarPost2CTL,ind,ind2)

        figure;hold on;
     VarPost2CTL(7)=83.1951/mean(mean(Experiment7.pitchACpostCTL(170:340,:))); % after removing outliers   
        figure;hold on;subplot(121);hold on;plot([1;2;4],[VarPre(ind);VarAPV(ind);VarPost2(ind)],'b-')
        plot([1;2;3;4],[VarPre(find(VarPost1>0));VarAPV(find(VarPost1>0));VarPost1(find(VarPost1>0));VarPost2(find(VarPost1>0))],'r-')
        plot(1,mean(VarPre(ind)),'k+','Markersize',25);plot(2,mean(VarAPV(ind)),'k+','Markersize',25);plot(4,mean(VarPost2(ind)),'k+','Markersize',25)
        plot(1,mean(VarPre(find(VarPost1>0))),'r+','Markersize',25);plot(2,mean(VarAPV(find(VarPost1>0))),'r+','Markersize',25);plot(3,mean(VarPost1(find(VarPost1>0))),'r+','Markersize',25);plot(4,mean(VarPost2(find(VarPost1>0))),'r+','Markersize',25)
        subplot(122);hold on;plot([1;2;4],[VarPreCTL(ind2);VarAPVCTL(ind2);VarPost2CTL(ind2)],'b-')  
        plot(1,mean(VarPreCTL(ind2)),'k+','Markersize',25);plot(2,mean(VarAPVCTL(ind2)),'k+','Markersize',25);plot(4,mean(VarPost2CTL(ind2)),'k+','Markersize',25)      
        plot([1;2;3;4],[VarPreCTL(find(VarPost1CTL>0));VarAPVCTL(find(VarPost1CTL>0));VarPost1CTL(find(VarPost1CTL>0));VarPost2CTL(find(VarPost1CTL>0))],'r-')   
        plot(1,mean(VarPreCTL(find(VarPost1CTL>0))),'r+','Markersize',25);plot(2,mean(VarAPVCTL(find(VarPost1CTL>0))),'r+','Markersize',25);plot(3,mean(VarPost1CTL(find(VarPost1CTL>0))),'r+','Markersize',25);plot(4,mean(VarPost2CTL(find(VarPost1CTL>0))),'r+','Markersize',25)