function p=NS248probability(deltaP,deltaT,T)
p=(0.5-deltaP)*[(1/(1+deltaT))*exp(-T/(1+deltaT))]+(0.5+deltaP)*[(1/(1-deltaT))*exp(-T/(1-deltaT))];
end
