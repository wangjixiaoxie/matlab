%purpose of this code is to go through sumbs and calculate
%the offz for each initrun
%for each control run
%and for each baseline runs
%save to a masterstruct suminitshs
%with initacz, initmuz, initoff, init_bsvl, 
%basacz,basmuz,basoff,bas_bsvl
%ctrlacz,ctrlmuz,ctrloff,ctrl_bsvl

figure
init.ac=[];
init.mu=[];
init.acer=[];
init.muer=[];

ctrl.ac=[];
ctrl.mu=[];
ctrl.acer=[];
ctrl.muer=[];

bas.ac=[];
bas.mu=[];
bas.acer=[];
bas.muer=[];



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
        
            
        acvls=acz(nt,initindcomb);
        acervls=aczer(nt,initindcomb);
        muvls=muz(nt, initindcomb);
        muervls=muzer(nt,initindcomb);
        
        if (ntind==1)
            init.ac=[init.ac acvls]
            init.mu=[init.mu muvls]
            init.acer=[init.acer  acervls]
            init.muer=[init.muer muervls]
        else
            ctrl.ac=[ctrl.ac acvls]
            ctrl.mu=[ctrl.mu muvls]
            ctrl.acer=[ctrl.acer  acervls]
            ctrl.muer=[ctrl.muer muervls]
        end
        
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
   
            bas.ac=[bas.ac acbasvls]
            bas.mu=[bas.mu mubasvls]
            bas.acer=[bas.acer  acbaservls]
            bas.muer=[bas.muer mubaservls]
   
        plot([acbasvls+acbaservls;acbasvls-acbaservls], [mubasvls;mubasvls],'Color','c');
        hold on;
        plot([acbasvls;acbasvls],[mubasvls+mubaservls;mubasvls-mubaservls],'Color','c');
   
   
   ylabel('muz');
   xlabel('acz');
end
   

bas.off=bas.ac-bas.mu;
ctrl.off=ctrl.ac-ctrl.mu;
%separate variables into init

initdown=[];
initup=[];
for ii=1:length(init.ac)
    if(init.ac(ii)>0)
        downtmp=init.ac(ii)-init.mu(ii);
        initdown=[initdown downtmp]
    else
        uptmp=init.ac(ii)-init.mu(ii);
        initup=[initup uptmp]
    end
end









    
