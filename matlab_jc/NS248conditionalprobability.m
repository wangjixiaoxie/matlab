function cp=NS248conditionalprobability(deltaP,deltaT,T)
y1=(0.5-deltaP);
y2=(0.5+deltaP);
x1=-1*trapz([(1/(1+deltaT))*exp(-T/(1+deltaT))]);
x2=-1*trapz([(1/(1-deltaT))*exp(-T/(1-deltaT))]);
cp=x1*y1+x2*y2;

