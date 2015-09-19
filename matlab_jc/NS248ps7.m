% PROBLEM 1
A.
H=-sum(p(x)*log2(p(x)))
% e.g. 2^-3=1/8 ...etc...
H=-(1/8*(-3)+1/2*(-1)+1/4*(-2)+1/8*(-3))
=1.75

B.
% Minimize the code by giving the most probable outcomes the shortest lengths
X=2 ---- 0 ----- 1*0.5
X=3 ---- 10 ---- 2*0.25
X=1 ---- 110 --- 3*0.125
X=4 ---- 1110 -- 4*0.125
% total = 1.875 --- slightly greater than H

% BUT we can do even better
X=2 ---- 0 ---- 1*0.5
X=3 ---- 10 --- 2*0.25
X=1 ---- 110 -- 3*0.125
X=2 ---- 111 -- 4*0.125
% total = 1.75 - same as H

C.
% Is it 2? -- 50% of the time it is 2, and you're done (0.5)*(1 question)
% Is it 3? -- 25% of the time it is 3, and you're done (0.25)*(2 questions)
% Is it 1? -- remaining 25% - you're done either way (0.25)*(3 questions)
    % because no means it has to be 4
Total = 1.75 = H

D. 
a. 
H  =  -(3*1/3*(log2(1/3)))  =  -log2(1/3)  =  1.5850
b.
X=1 --- 0
X=2 --- 10
X=3 --- 11
total  =  (1+2+2)/3  =  5/3  =  1.6667>1.5850
c.
% Is it 1? -- (1/3)*(1 question)
% Is it 2? -- (2/3)*(2 questions)
Total = 5/3 = 1.6667>1.5850


% PROBLEM 2  
% I AM UNSURE ABOUT THE T values for numerical integration
a.
H=-INTEGRAL(p(x)*log2(p(x))
% Spike at time T can happen either due to a weak flash or a strong flash
function p=NS248probability(deltaP,deltaT,T)
p=(0.5-deltaP)*[(1/(1+deltaT))*exp(-T/(1+deltaT))]+(0.5+deltaP)*[(1/(1-deltaT))*exp(-T/(1-deltaT))];
end

deltaTvals=[0.01:0.01:0.99];
deltaPvals=[0.01:0.01:0.49];
for i=1:length(deltaTvals)
    deltaT=deltaTvals(i);
    for j=1:length(deltaPvals)
        deltaP=deltaPvals(j);
        T=[0:0.1:100];  % T=100 gives the exponential enough time to decay to zero (same answer for T=1000);
                        % KEEP THE T-step value constant
        p=NS248probability(deltaP,deltaT,T);
        entropy(i,j)=-trapz(p.*log2(p));
    end
end
% Entropy is greatest when deltaT is close to zero and decreases as deltaP approaches 0.5
figure;imagesc(entropy(6:89,:))
% PS7figure1CORRECTVERSION.png --- in /cardinal/NS248 folder

B.
E=I/H(X)=[H(X)-H(X|Y)]/H(X)
% we know H(X) from the pvs section
H(X|Y)=sum[p(stimulusN)*(-spike distn given stimulusN*log2(spike distn given stimulusN)]
    stimulusN has two values ---> weak flash and strong flash
function cp=NS248conditionalprobability(deltaP,deltaT,T)
y1=(0.5-deltaP);
y2=(0.5+deltaP);
x1=-1*trapz([(1/(1+deltaT))*exp(-T/(1+deltaT))]);
x2=-1*trapz([(1/(1-deltaT))*exp(-T/(1-deltaT))]);
cp=x1*y1+x2*y2;
end
    
deltaTvals=[0.01:0.01:0.99];
deltaPvals=[0.01:0.01:0.49];
for i=1:length(deltaTvals)
    deltaT=deltaTvals(i);
    for j=1:length(deltaPvals)
        deltaP=deltaPvals(j);
        T=[0:0.1:100];  % T=100 gives the exponential enough time to decay to zero (same answer for T=1000);
        entro=entropy(i,j);
        cp=NS248conditionalprobability(deltaP,deltaT,T);
        E(i,j)=(entro-cp)/entro;
    end
end
figure;imagesc(E)
% PS7figure2.png

C.
pspikeA=1-deltaP;
pspikeB=0.5+deltaP;
entropy2=-pspikeA*log2(pspikeA)-pspikeB*log2(pspikeB);
cp2=