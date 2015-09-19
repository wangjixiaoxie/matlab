function doitnow(p716lid,f716lid,p716B,f716B,avgnote)
CVlido1=jc_plotabunch725(p716lid(1:250),f716lid(1:250),900,avgnote,400,600);
CVacsf1=jc_plotabunch725(p716B(1:250),f716B(1:250),900,avgnote,400,600);
for i=0:4
    g=CVlido1(i*50+1:(i+1)*50);
    h=CVacsf1(i*50+1:(i+1)*50);
    CVL(i+1)=std(g)/mean(g);
    CVA(i+1)=std(h)/mean(h);
end
figure; hold on; plot(CVL); plot(CVA,'r')
    