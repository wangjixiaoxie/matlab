function [xvals,propvals]=PSDcontribution(Data,Datakeep,index,width)
% e.g. PSDcontribution(BFLesion,[1 2 3 4 9],[950 450 550 450 300],25);
% Data - should be formatted like BFLesion
% index - middle of note
% width - width of chunk of the note in milliseconds



ptstart=index-(width*8)/2+1;
ptend=index+(width*8)/2;

%%%%% METHOD 1: Relative PSD at different frequencies
% Calculate the power spectral density across your window
    for k=1:length(Datakeep)
        i=Datakeep(k);
        [x,psdsPRE(k).data]=jcpsd2(Data(i).UDpre20_1024(ptstart(k):ptend(k),:),8000);
        [x,psdsPOST(k).data]=jcpsd2(Data(i).UDpost20_1024(ptstart(k):ptend(k),:),8000);
    end
    ps=find(x>5);
% Calculate the proportion of reduction in power -> (pre-post)/pre
    for j=1:length(psdsPRE)
        mvPRE(j,ps)=median(psdsPRE(j).data(:,ps));
        mvPOST(j,ps)=median(psdsPOST(j).data(:,ps));
    end
    xvals=x(ps);
    propvals=(mvPRE-mvPOST)./mvPRE;
    
% figure;hold on
% for i=1:length(psdsPREBF)
%     plot(xB,median(psdsPOSTBF(i).data)./median(psdsPREBF(i).data),'r')
% end
% for i=1:length(psdsPREBF)
%     plot(xB',median(psdsPREBF(i).data),'*')
%     plot(xB',median(psdsPOSTBF(i).data),'*')
% end
% figure;hold on
% for i=1:length(psdsPREBF)
%     plot(xB10',median(psdsPREBF10ms(i).data),'k')
%     plot(xB10',median(psdsPOSTBF10ms(i).data),'r')
% end
% for i=1:length(psdsPREBF)
%     plot(xB10',median(psdsPREBF10ms(i).data),'*')
%     plot(xB10',median(psdsPOSTBF10ms(i).data),'*')
% end
% 
% figure;hold on
% for i=1:length(psdsPRE)
% plot(xZ',median(psdsPRE(i).data)./median(psdsPOST(i).data),'r')
% end
% 
% clear acorr
% for i=1:length(successes)
%     clear autoc
%     for j=1:size(BFLesion(i).UDpost10_512,2)
%         autoc(j,:)=xcorr(BFLesion(i).UDpost10_512(:,j));
%     end
%     acorr(i,:)=median(autoc);
% end





