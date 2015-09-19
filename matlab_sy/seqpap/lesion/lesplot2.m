%include evren's data:

%plot lesion data
% lesplot(bs,bnum)
%PLOT1 (321) is pitch pre post learning, 
%PLOT2 (322)  is sequence pre post learning
%PLOT3 (323)  pitch learning vs. sequence learning for matched pairs
%PLOT4 (325:6)  is bar graph comparing learning - with individual
%comparisons...


%third plot is pitch learning vs. sequence learning for paired birds
%5th and 6th plots are bar graphs of net learning
function [outmt]=lesplot2()


%my paired lesion data:
cd /oriole7/dir1/lesionsumdata/
load sumdata.mat
npair=length(outmt);

%evren's data:
cd /oriole5/evrenlesiondata/sumdata
load postshiftout.mat
mrklst={'s' '+' 'x','h','*'}
colors={'r' 'k' 'c' 'b' 'm'}

combmt.ptpre=zeros(2, npair);
combmt.ptpst=zeros(2,npair);
combmt.seqpre=zeros(2,npair);
combmt.seqpst=zeros(2,npair);

figure
subplot(321)

%first combine data for my data
for ii=1:length(outmt)
%first plot pitch pre
    %first row is baspitch
    combmt.ptpre(1,ii)= outmt(ii).prelespt.freqbas
    %second row is wn pitch
    combmt.ptpre(2,ii)=outmt(ii).prelespt.freqsh
    if(length(outmt(ii).pstlespt)>1)
        pstpt=outmt(ii).pstlespt(1)
    else
        pstpt=outmt(ii).pstlespt
    end
    combmt.ptpst(1,ii)= pstpt.freqbas
    combmt.ptpst(2,ii)= pstpt.freqsh
  
 
end

%combine evren's data:


%plot both data in black for pre and red color
for ii=1:npair
    plot(combmt.ptpre(1,ii),combmt.ptpst(1,ii),'Marker',mrklst{ii},'Color','k')
    hold on;
    plot(combmt.ptpre(2,ii),combmt.ptpst(2,ii),'Marker',mrklst{ii},'Color','r')
end
    
%     vls2=outmt(ii).prelespt.freqsh-outmt(ii).prelespt.freqbas
% %     vls3=outmt(ii).prelespt.freqrec-outmt(ii).prelespt.freqbas
%     if(outmt(ii).drxn=='dn')
%             vls2=-vls2
% %             vls3=-vls3
%     end
%     
%     plot([1 2 ],[vls1 vls2 ],colors{ii})
%     hold on;
%     plot([1 2 ],[vls1 vls2 ],'o','Color',colors{ii})
%     combmt.ptpre=[combmt.ptpre vls2-vls1];
%      
%     
% ax(1)=subplot(232)
%     if(~isempty(outmt(ii).pstlespt))
%         vls1=0
%        
        
%             vls2=crpt.freqsh-crpt.freqbas
% %             vls3=crpt.freqrec-crpt.freqbas
%         if(outmt(ii).drxn=='dn')
%             vls2=-vls2
% %             vls3=-vls3
%         end
%         
%         plot([1 2 ],[vls1 vls2 ],colors{ii})
%         hold on;
%         plot([1 2 ],[vls1 vls2 ],'o','Color',colors{ii})
%         combmt.ptpst=[combmt.ptpst vls2-vls1];
%         
%     end

 
 
%     ax(2)=subplot(234)
% 
%     vls=[outmt(ii).preles_sq.seqbas outmt(ii).preles_sq.seqsh ]
% %      vls=[outmt(ii).preles_sq.seqbas outmt(ii).preles_sq.seqsh outmt(ii).preles_sq.seqrec]
%     plot([1 2 ],1-vls,colors{ii})
%     hold on;
%     plot([1 2 ],1-vls,'o','Color',colors{ii})
%      combmt.seqpre=[combmt.seqpre vls(1)-vls(2)];
%         
%     
% ax(2)=subplot(235)
% %      vls=[outmt(ii).pstles_sq.seqbas(end) outmt(ii).pstles_sq.seqsh(end) outmt(ii).pstles_sq.seqrec(end)]
%       vls=[outmt(ii).pstles_sq.seqbas(end) outmt(ii).pstles_sq.seqsh(end) ]
%      plot([1 2],1-vls,colors{ii})
%      hold on;
%      plot([1 2 ],1-vls,'o','Color',colors{ii})
%      combmt.seqpst=[combmt.seqpst vls(1)-vls(2)]
% end
% axtt=subplot(233)
%     mn(1)=mean(combmt.ptpre);
%     er(1)=std(combmt.ptpre)./sqrt(length(combmt.ptpre));
%     mn(2)=mean(combmt.ptpst);
%     er(2)=std(combmt.ptpst)./sqrt(length(combmt.ptpst));
%     makebar(combmt.ptpre,combmt.ptpst,mn,er);
% axtt=subplot(236)
%     mn(1)=mean(combmt.seqpre);
%     er(1)=std(combmt.seqpre)./sqrt(length(combmt.seqpre));
%     mn(2)=mean(combmt.seqpst);
%     er(2)=std(combmt.seqpst)./sqrt(length(combmt.seqpst));
%     makebar(combmt.seqpre,combmt.seqpst,mn,er);
    
    
    function []=makebar(prevl,pstvls,mn,er)
        plot([ones(1,4);2*ones(1,4)],[prevl ;pstvls],'k')
        hold on;
        plot(ones(1,4),prevl,'ko');
        plot(2*ones(1,4),pstvls,'ko')
        hold on;
        bar([1 2],[mn(1) mn(2)],'FaceColor','none');
        hold on;
        plot([1  1],[mn(1)+er(1) mn(1)-er(1)])
        plot([2  2],[mn(2)+er(2) mn(2)-er(2)])   
