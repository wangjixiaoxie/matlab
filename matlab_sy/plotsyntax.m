%general plotting program for syntax transitions.  In progress.

days=[1:12]

for ii=1:length(days)
    dv=days(ii);
    efrac(ii)=synstruct{dv}.e/synstruct{dv}.ntrans;
    afrac(ii)=synstruct{dv}.a/synstruct{dv}.ntrans;
    ifrac(ii)=synstruct{dv}.i/synstruct{dv}.ntrans;
    cfrac(ii)=synstruct{dv}.c/synstruct{dv}.ntrans;
    ofrac(ii)=1-efrac(ii)-afrac(ii)-ifrac(ii)-cfrac(ii);
    day(ii)=synstruct{dv}.day
end

ofrac=ofrac+efrac
figure
x1=floor(switchdts(1))-day(1)-.9
x2=floor(switchdts(2))-day(2)+.1
y1=0
y2=1;
%ax=fill([x1 x2 x2 x1],[y1 y1 y2 y2],[1 .806 .817]);
hold on;
%plot(day(2:11)-day(2), efrac(2:11),'k-o','Linewidth',2);hold on;
plot(day(2:11)-day(2),afrac(2:11),'b-o')
plot(day(2:11)-day(2),ofrac(2:11),'k-o')
plot(day(2:11)-day(2),ifrac(2:11),'g-o')
plot(day(2:11)-day(2),cfrac(2:11),'r-o')

axis([0 14 0 1]);