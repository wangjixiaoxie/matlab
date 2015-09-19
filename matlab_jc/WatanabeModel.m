
figure;hold on;
for j=1:12
integC=0.005*5*j;
for i=1:1000
a=1; % produce repeated syllable
b=0; % produce subsequent syllable
c=1;
count=1;
d=1;
while d>0
    a(count)=1; % fires first, no memory
    c(count+1)=c(count)+a(count)*integC*rand; % fires second, has memory
    b(count+1)=b(count)+c(count+1)-a(count); % receives input third, no memory
    a(count+1)=a(count)-b(count+1);
    count=count+1;
    d=a(count)-b(count);
end
lengthdist(i)=count-1;
end
[b,a]=hist(lengthdist,[1:1:20]);
subplot(4,3,j);stairs(b/length(lengthdist))
xlim([0 20])
ylim([0 0.65])
end