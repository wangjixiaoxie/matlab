

pu24ccDIR(1)=1;
for i=8:10
    a=corrcoef(mean(ZFLesionInternote(7).pitchDIR(ZFLesionInternote(7).cntr-80:ZFLesionInternote(7).cntr+80,:)),mean(ZFLesionInternote(i).pitchDIR(ZFLesionInternote(i).cntr-80:ZFLesionInternote(i).cntr+80,:)));
    pu24ccDIR(i-6)=a(2);
end
    a=corrcoef(mean(ZFLesionInternote(8).pitchDIR(ZFLesionInternote(8).cntr-80:ZFLesionInternote(8).cntr+80,:)),mean(ZFLesionInternote(9).pitchDIR(ZFLesionInternote(9).cntr-80:ZFLesionInternote(9).cntr+80,:)));
    pu24ccDIR(5)=a(2);
    a=corrcoef(mean(ZFLesionInternote(8).pitchDIR(ZFLesionInternote(8).cntr-80:ZFLesionInternote(8).cntr+80,:)),mean(ZFLesionInternote(10).pitchDIR(ZFLesionInternote(10).cntr-80:ZFLesionInternote(10).cntr+80,:)));
    pu24ccDIR(6)=a(2);
    a=corrcoef(mean(ZFLesionInternote(9).pitchDIR(ZFLesionInternote(9).cntr-80:ZFLesionInternote(9).cntr+80,:)),mean(ZFLesionInternote(10).pitchDIR(ZFLesionInternote(10).cntr-80:ZFLesionInternote(10).cntr+80,:)));
    pu24ccDIR(7)=a(2);
    
 figure;plot(pu24ccDIR,'g')
hold on;plot(pu24ccPRE,'b')
hold on;plot(pu24ccDIR,'r')
hold on;plot(pu24ccDIR,'g')
hold on;plot(pu24ccDIR,'*','Color','g')
hold on;plot(pu24ccDIR,'*','Color','r')
hold on;plot(pu24ccPRE,'*','Color','b')

    

corrcoef(mean(RepeatsBF(3).pitchs2a(400:500,:)),mean(RepeatsBF(3).pitchs2b(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2a(400:500,:)),mean(RepeatsBF(3).pitchs2c(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2a(400:500,:)),mean(RepeatsBF(3).pitchs2d(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2a(400:500,:)),mean(RepeatsBF(3).pitchs2e(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2d(400:500,:)),mean(RepeatsBF(3).pitchs2e(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2c(400:500,:)),mean(RepeatsBF(3).pitchs2e(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2c(400:500,:)),mean(RepeatsBF(3).pitchs2d(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2b(400:500,:)),mean(RepeatsBF(3).pitchs2d(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2b(400:500,:)),mean(RepeatsBF(3).pitchs2c(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2b(400:500,:)),mean(RepeatsBF(3).pitchs2e(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2aD(400:500,:)),mean(RepeatsBF(3).pitchs2bD(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2aD(400:500,:)),mean(RepeatsBF(3).pitchs2cD(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2aD(400:500,:)),mean(RepeatsBF(3).pitchs2dD(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2aD(400:500,:)),mean(RepeatsBF(3).pitchs2eD(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2bD(400:500,:)),mean(RepeatsBF(3).pitchs2cD(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2bD(400:500,:)),mean(RepeatsBF(3).pitchs2dD(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2bD(400:500,:)),mean(RepeatsBF(3).pitchs2eD(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2dD(400:500,:)),mean(RepeatsBF(3).pitchs2eD(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2cD(400:500,:)),mean(RepeatsBF(3).pitchs2eD(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2cD(400:500,:)),mean(RepeatsBF(3).pitchs2dD(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2cINA(400:500,:)),mean(RepeatsBF(3).pitchs2dINA(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2cINA(400:500,:)),mean(RepeatsBF(3).pitchs2eINA(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2dINA(400:500,:)),mean(RepeatsBF(3).pitchs2eINA(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2bINA(400:500,:)),mean(RepeatsBF(3).pitchs2eINA(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2bINA(400:500,:)),mean(RepeatsBF(3).pitchs2dINA(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2bINA(400:500,:)),mean(RepeatsBF(3).pitchs2cINA(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2aINA(400:500,:)),mean(RepeatsBF(3).pitchs2cINA(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2aINA(400:500,:)),mean(RepeatsBF(3).pitchs2bINA(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2aINA(400:500,:)),mean(RepeatsBF(3).pitchs2dINA(400:500,:)))
corrcoef(mean(RepeatsBF(3).pitchs2aINA(400:500,:)),mean(RepeatsBF(3).pitchs2eINA(400:500,:)))
    


figure;hold on
for i=1:5
    plot([1 RepeatsBF(i).ccPRE],'b')
    plot([1 RepeatsBF(i).ccPRE],'*','Color','b')
        plot([1 RepeatsBF(i).ccPOST],'r')
    plot([1 RepeatsBF(i).ccPOST],'*','Color','r')
        plot([1 RepeatsBF(i).ccDIR],'g')
    plot([1 RepeatsBF(i).ccDIR],'*','Color','g')
end
for i=3
    plot([1 RepeatsBF(i).ccPRE2],'b')
    plot([1 RepeatsBF(i).ccPRE2],'*','Color','b')
        plot([1 RepeatsBF(i).ccPOST2],'r')
    plot([1 RepeatsBF(i).ccPOST2],'*','Color','r')
        plot([1 RepeatsBF(i).ccDIR2],'g')
    plot([1 RepeatsBF(i).ccDIR2],'*','Color','g')
end
