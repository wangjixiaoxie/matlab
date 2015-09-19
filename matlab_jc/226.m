[vals313,trigs]=triglabel('batch225files','a',1,1,0,1);
fvs313=findwnoteJC('batch225notes','a','','',0,[3000 4200],8000,0,'obs0',0);
tvals313=timing(vals225,fvs225);
FFvals313=evtaf_freq('batch225files',[3200 4000],'a',128,'obs0',1,1);
FFvals313=FFvals313(:,2);
figure;plot([tvalsBaseline tvalsTRANS tvalsAPV tvalsPOST,tvalsPOST2,tvalsPROBE,tvalsTRANS2,tvalsAPV2,tvalsAPV5,tvalsPOST3,tvalsPOST4],[FFvalsBaseline;FFvalsTRANS;FFvalsAPV;FFvalsPOST;FFvalsPOST2;FFvalsPROBE;FFvalsTRANS2;FFvalsAPV2;FFvalsAPV5;FFvalsPOST3;FFvalsPOST4],'*','Color','r')
hold on;plot([tvalsBaseline tvalsTRANS tvalsAPV tvalsPOST,tvalsPOST2,tvalsPROBE,tvalsTRANS2,tvalsAPV2,tvalsAPV5],[FFvalsBaseline;FFvalsTRANS;FFvalsAPV;FFvalsPOST;FFvalsPOST2;FFvalsPROBE;FFvalsTRANS2;FFvalsAPV2;FFvalsAPV5],'*','Color','k')
hold on;plot([tvalsBaseline tvalsTRANS tvalsAPV tvalsPOST,tvalsPOST2,tvalsPROBE,tvalsTRANS2,tvalsAPV2],[FFvalsBaseline;FFvalsTRANS;FFvalsAPV;FFvalsPOST;FFvalsPOST2;FFvalsPROBE;FFvalsTRANS2;FFvalsAPV2],'*','Color','k')
hold on;plot([tvalsBaseline tvalsTRANS tvalsAPV tvalsPOST,tvalsPOST2,tvalsPROBE,tvalsTRANS2],[FFvalsBaseline;FFvalsTRANS;FFvalsAPV;FFvalsPOST;FFvalsPOST2;FFvalsPROBE;FFvalsTRANS2],'*','Color','g')
hold on;plot([tvalsBaseline tvalsTRANS tvalsAPV tvalsPOST,tvalsPOST2,tvalsPROBE],[FFvalsBaseline;FFvalsTRANS;FFvalsAPV;FFvalsPOST;FFvalsPOST2;FFvalsPROBE],'*','Color','b')
hold on;plot([tvalsBaseline tvalsTRANS tvalsAPV tvalsPOST,tvalsPOST2],[FFvalsBaseline;FFvalsTRANS;FFvalsAPV;FFvalsPOST;FFvalsPOST2],'*','Color','r')
hold on;plot([tvalsBaseline tvalsTRANS tvalsAPV],[FFvalsBaseline;FFvalsTRANS;FFvalsAPV],'*','Color','k')
hold on;plot([tvalsBaseline tvalsTRANS],[FFvalsBaseline;FFvalsTRANS],'*','Color','g')
hold on;plot([tvalsBaseline],[FFvalsBaseline],'*')