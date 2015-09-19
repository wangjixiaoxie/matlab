[valsREC,trigs]=triglabel('batchRECfiles','a',1,1,0,1);
fvsREC=findwnoteJC('batchRECnotes','a','','',0,[3000 4200],8000,0,'obs0',0);
tvalsREC=timing(valsREC,fvsREC);
FFvalsREC=evtaf_freq('batchRECfiles',[3200 4000],'a',128,'obs0',1,1);
FFvalsREC=FFvalsREC(:,2);
figure;plot([tvalsBaseline tvalsTRANS tvalsAPV tvalsPOST,tvalsPOST2,tvalsPROBE,tvalsTRANS2,tvalsAPV2,tvalsAPV5,tvalsPOST3,tvalsPOST4,tvals313,tvalsREC],[FFvalsBaseline;FFvalsTRANS;FFvalsAPV;FFvalsPOST;FFvalsPOST2;FFvalsPROBE;FFvalsTRANS2;FFvalsAPV2;FFvalsAPV5;FFvalsPOST3;FFvalsPOST4;FFvals313;FFvalsREC],'*','Color','r')
hold on;plot([tvalsBaseline tvalsTRANS tvalsAPV tvalsPOST,tvalsPOST2,tvalsPROBE,tvalsTRANS2,tvalsAPV2,tvalsAPV5,tvalsPOST3,tvalsPOST4,tvals313],[FFvalsBaseline;FFvalsTRANS;FFvalsAPV;FFvalsPOST;FFvalsPOST2;FFvalsPROBE;FFvalsTRANS2;FFvalsAPV2;FFvalsAPV5;FFvalsPOST3;FFvalsPOST4;FFvals313],'*','Color','k')
hold on;plot([tvalsBaseline tvalsTRANS tvalsAPV tvalsPOST,tvalsPOST2,tvalsPROBE,tvalsTRANS2,tvalsAPV2,tvalsAPV5,tvalsPOST3,tvalsPOST4],[FFvalsBaseline;FFvalsTRANS;FFvalsAPV;FFvalsPOST;FFvalsPOST2;FFvalsPROBE;FFvalsTRANS2;FFvalsAPV2;FFvalsAPV5;FFvalsPOST3;FFvalsPOST4],'*','Color','r')
hold on;plot([tvalsBaseline tvalsTRANS tvalsAPV tvalsPOST,tvalsPOST2,tvalsPROBE,tvalsTRANS2,tvalsAPV2,tvalsAPV5],[FFvalsBaseline;FFvalsTRANS;FFvalsAPV;FFvalsPOST;FFvalsPOST2;FFvalsPROBE;FFvalsTRANS2;FFvalsAPV2;FFvalsAPV5],'*','Color','k')
hold on;plot([tvalsBaseline tvalsTRANS tvalsAPV tvalsPOST,tvalsPOST2,tvalsPROBE,tvalsTRANS2,tvalsAPV2],[FFvalsBaseline;FFvalsTRANS;FFvalsAPV;FFvalsPOST;FFvalsPOST2;FFvalsPROBE;FFvalsTRANS2;FFvalsAPV2],'*','Color','k')
hold on;plot([tvalsBaseline tvalsTRANS tvalsAPV tvalsPOST,tvalsPOST2,tvalsPROBE,tvalsTRANS2],[FFvalsBaseline;FFvalsTRANS;FFvalsAPV;FFvalsPOST;FFvalsPOST2;FFvalsPROBE;FFvalsTRANS2],'*','Color','g')
hold on;plot([tvalsBaseline tvalsTRANS tvalsAPV tvalsPOST,tvalsPOST2,tvalsPROBE],[FFvalsBaseline;FFvalsTRANS;FFvalsAPV;FFvalsPOST;FFvalsPOST2;FFvalsPROBE],'*','Color','b')
hold on;plot([tvalsBaseline tvalsTRANS tvalsAPV tvalsPOST,tvalsPOST2],[FFvalsBaseline;FFvalsTRANS;FFvalsAPV;FFvalsPOST;FFvalsPOST2],'*','Color','r')
hold on;plot([tvalsBaseline tvalsTRANS tvalsAPV],[FFvalsBaseline;FFvalsTRANS;FFvalsAPV],'*','Color','k')
hold on;plot([tvalsBaseline tvalsTRANS],[FFvalsBaseline;FFvalsTRANS],'*','Color','g')
hold on;plot([tvalsBaseline],[FFvalsBaseline],'*')

% vals313=control shift
% valsREC=recovery from control shift