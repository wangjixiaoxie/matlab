% minstim=.062
% maxstim=.004
clear stcr
muinds=[51]
for ii=1:length(muinds)
    run_num=muinds(ii);
    bt=avls.cvl{run_num}
    pathvl=avls.pvls{run_num}
    if (isfield(avls,'baspath'))
        cmd=['cd ' avls.baspath pathvl]
    else
        cmd=['cd ' pathvl]
    end

    eval(cmd);
    
    cmd=['load ' bt '.mat']
    eval(cmd);
    
    
    
    clear minstim
    clear maxstim
    clear diffcontour

    delay=avls.del{run_num}
    
    tmbuff=0.1
    contouroffset=avls.con_tbinshft;
    contourtms=[33 363]
    mncontour=mean(contours(contourtms(1):contourtms(2),crctind),2);

    for jj=1:length(contours(1,:))
        diffcontour{jj}=contours(contourtms(1):contourtms(2),jj)-mncontour;
        
    end

%restrict the size of the contour from .07 to .12

    tmdiff=pitchtms(2)-pitchtms(1);

    initmat=zeros(1000,length(diffcontour));


     for jj=1:length(fvst)
       stcr(jj)=fvst(jj).STIMTIME;
            if(stcr(jj)~=-1)
                
          
                tmshft(jj)=(stcr(jj)/1000)+delay;
                adjtmshft(jj)=pitchtms(contourtms(1))-tmshft(jj);
                ptshft(jj)=round(adjtmshft(jj)/tmdiff);
                
                
%                 ptshft can be negative
             
                %set up an extra buffer, 20ms, and start times at a negative
                %value.
                ptbuff=tmbuff/tmdiff;
              
                initmat(floor(ptshft(jj)+ptbuff):floor(ptbuff+ptshft(jj)+length(diffcontour{jj})-1),jj)=diffcontour{jj}(1:end);
            
            end
   end
   mndiff=mean(initmat(:,crfbind),2)-mean(initmat(:,crctind),2);
   mat_tms=-tmbuff:tmdiff:tmdiff*length(initmat(:,1))-tmbuff;
   mat_tms=mat_tms(1:end-1);
   %generate mean and standard deviation ind by ind)
   for indvl=1:length(initmat(:,1))
       nozeroind=find(initmat(indvl,crfbind)~=0);
       if(~isempty(nozeroind))
            if(length(nozeroind)>10)
                mnout(indvl)=mean(initmat(indvl,crfbind(nozeroind)));
                ster(indvl)=std(initmat(indvl,crfbind(nozeroind))./sqrt(length(nozeroind)));
            else   
                mnout(indvl)=0;
                ster(indvl)=0;
            end
       else
         mnout(indvl)=0;
         ster(indvl)=0;
       end
   end
         
       
   
   
%    cmd=['save -append ' bt '.mat mat_tms initmat stcr diffcontour mnout ster']  ;
%     eval(cmd);
   
end
   