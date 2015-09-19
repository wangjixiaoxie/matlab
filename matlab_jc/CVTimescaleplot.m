function CVTimescaleplot(Alldata2)
 %%%%%%%%%%%
 % CVchange vs. TimescaleChange plot
 %%%%%%%%%%
 
 g=[1 4 5 6 10 11 15 16];
for i=1:length(g)
    n=g(i);
    midd=(Alldata2(n).basevaroffset-Alldata2(n).basevaronset)./2;
    AC=std(Alldata2(n).baselineAC(Alldata2(n).basevaronset:Alldata2(n).basevaroffset,:)')./mean(Alldata2(n).baselineAC(Alldata2(n).basevaronset:Alldata2(n).basevaroffset,:)');
    INA=std(Alldata2(n).baselineINA(Alldata2(n).basevaronset:Alldata2(n).basevaroffset,:)')./mean(Alldata2(n).baselineINA(Alldata2(n).basevaronset:Alldata2(n).basevaroffset,:)');
    CVac(i)=mean(AC(midd-40:midd+40));
    CVina(i)=mean(INA(midd-40:midd+40));
end
[slopesAC,slopesINA]=xcorranalyINA(Alldata2);
slopesratio=slopesAC./slopesINA;
g=[1 4 5 6 10 11 15 16];
for i=1:length(g)
n=g(i);
midd=(Alldata2(n).basevaroffset-Alldata2(n).basevaronset)./2;
AC=std(Alldata2(n).baselineAC(Alldata2(n).basevaronset:Alldata2(n).basevaroffset,:)')./mean(Alldata2(n).baselineAC(Alldata2(n).basevaronset:Alldata2(n).basevaroffset,:)');
INA=std(Alldata2(n).baselineINA(Alldata2(n).basevaronset:Alldata2(n).basevaroffset,:)')./mean(Alldata2(n).baselineINA(Alldata2(n).basevaronset:Alldata2(n).basevaroffset,:)');
CVac(i)=mean(AC(midd-40:midd+40));
CVina(i)=mean(INA(midd-40:midd+40));
end
CVratio=CVac./CVina;
CVratio=CVratio(1:8);
figure;plot(1./slopesratio,CVratio,'*')
for i=1:200
x(i)=0+0.01*i;
y(i)=+0.01*i;
end
hold on;plot(x,y)
ksim=polyfit((1./slopesratio),CVratio,1)
0.254*2.4
ans+0.8770
for i=1:200
xx(i)=i*0.01;
yy(i)=xx(i)*0.254+0.877;
end
