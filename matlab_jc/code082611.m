for i=1:100000
% Model 1 - amplification with a bit of jitter
HVC=(rand(1,11)*11)-5.5;
LMAN=HVC*rand*3;
RA1(i)=sum(HVC+LMAN);
% Model 1B - learning
HVC=(rand(1,11)*11)-5.5;
LMAN=abs(HVC*4*rand);
RA1B(i)=sum(HVC+LMAN);
% Model 2 - uncorrelated influences
HVC=(rand(1,11)*11)-5.5;
LMAN=(rand(1,11)*11)-5.5;
RA2(i)=sum(HVC+4*LMAN);
% Model 2B - learning
HVC=(rand(1,11)*11)-5.5;
LMAN=abs((rand(1,11)*11)-5.5);
RA2B(i)=sum(HVC+4*LMAN);
end
figure;hist(RA1B,[-150:10:250])
figure;hist(RA2B,[-150:10:150])

% 
round(rand(1,11)) % random activation of 11 neurons

 