%plotinitiind.m
%switch to STANDARD DEVIATION
figure
for ii=1:length(sumbs)
   ii
    acz=sumbs(ii).acz
   muz=sumbs(ii).muz
   aczer=sumbs(ii).acerrz
   muzer=sumbs(ii).muerrz
   initindcomb=[];
   for jj=1:length(sumbs(ii).initind)
       if(~isempty(sumbs(ii).initind{jj}))
            initindcomb=[initindcomb ;sumbs(ii).initind{jj}(1)]
       end
       basind=[sumbs(ii).basruns]
    end
        %switch this to all inds, not just targetind
   allnote=sumbs(ii).allnote;
%    plot(acz(ntind, initindcomb), muz(ntind,initindcomb),'o');
%    hold on;
   for ntind=1:length(allnote)
       nt=allnote(ntind);
        if ntind==1
           col='k' 
        else
            col='r'
        end
            
        acvls=acz(nt,initindcomb);
        acervls=aczer(nt,initindcomb);
        muvls=muz(nt, initindcomb);
        muervls=muzer(nt,initindcomb);
        
        
        
        plot([acvls+acervls;acvls-acervls], [muvls;muvls],'Color',col);
        hold on;
        plot([acvls;acvls],[muvls+muervls;muvls-muervls],'Color',col);
        plot([-10 10],[-10 10],'k-')
   
   
   end
   
   %ADD BASELINE
   nt=allnote(1);
   acbasvls=acz(nt,basind);
   acbaservls=aczer(nt,basind);
   mubasvls=muz(nt,basind);
   mubaservls=muzer(nt,basind);
   
   
        plot([acbasvls+acbaservls;acbasvls-acbaservls], [mubasvls;mubasvls],'Color','c');
        hold on;
        plot([acbasvls;acbasvls],[mubasvls+mubaservls;mubasvls-mubaservls],'Color','c');
   
   
   ylabel('muz');
   xlabel('acz');
end
   
    
