function plotPreSleepFigure(Permanent,VarAPV,VarAPVCTL,indDN,Learn1,Learn2,mpitchpreTarg,ind3,Acute,FFpermanent1,VarPreCTL,VarPost1CTL,VarPost1,VarPost2CTL,VarPost2,VarPre)
subplot(131);hold on; % Learning
plot(Learn1(indDN)./mpitchpreTarg(indDN),'b')
plot(Learn1(indDN)./mpitchpreTarg(indDN),'b.')
plot(Learn2(indDN)./mpitchpreTarg(indDN),'r')
plot(Learn2(indDN)./mpitchpreTarg(indDN),'r.')
plot([0.5 9.5],[0 0],'k')
xlim([0.5 8.5])
subplot(132);hold on; % Var, control notes
plot(VarPost1CTL(ind3)./VarPreCTL(ind3),'b')
plot(VarPost1(indDN)./VarPre(indDN),'k')
plot(VarAPV(indDN)./VarPre(indDN),'r')
plot(VarAPVCTL(indDN)./VarPreCTL(indDN),'g')
plot(VarPost2CTL(ind3)./VarPreCTL(ind3),'b')
plot(VarPost2(indDN)./VarPre(indDN),'k')
plot(VarPost1CTL(ind3)./VarPreCTL(ind3),'b.')
plot(VarPost1(indDN)./VarPre(indDN),'k.')
plot(VarPost2CTL(ind3)./VarPreCTL(ind3),'b.')
plot(VarPost2(indDN)./VarPre(indDN),'k.')
plot(VarAPV(indDN)./VarPre(indDN),'r.')
plot(VarAPVCTL(indDN)./VarPreCTL(indDN),'g.')
plot([0.5 9.5],[1 1],'k')
xlim([0.5 8.5])
subplot(133);hold on; % Acute
plot(Acute(ind3)./mpitchpreTarg(ind3),'r')
plot(Acute(ind3)./mpitchpreTarg(ind3),'r.')
plot(FFpermanent1(ind3)./mpitchpreTarg(ind3),'b')
plot(FFpermanent1(ind3)./mpitchpreTarg(ind3),'b.')
plot(Permanent(ind3)./mpitchpreTarg(ind3),'b')
plot(Permanent(ind3)./mpitchpreTarg(ind3),'k.')
plot([0.5 6.5],[0 0],'k')
xlim([0.5 6.5])
