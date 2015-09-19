%this is a general interpolation function which takes a series of x
   %values and a series of y values and outputs y values for any gaps in
   %the x values 
   %for example x=1 2 4
   %            y=8 8 10
   %output is 8 8 9 10
   
   function [intout]=interp(yvls,xvls);             
       dys=xvls;
      
       intout(1)=yvls(1);
       for dyind=2:length(dys)
            dyvlcur=dys(dyind);
            intout(dyvlcur)=combeff(dyind);
            dyvlpre=yvls(dyind-1);
            diffvl=dyvlcur-dyvlpre;

            combeffpst=yvls(dyind)
            combeffpre=yvls(dyind-1);
            %diffdys 1 no interpr
            if(diffvl>1)
                for diffind=1:(diffvl-1);
                    prewt=((diffvl-diffind)/diffvl);
                    postwt=diffind/diffvl;
                    intout(dyvlpre+diffind)=prewt*combeffpre+postwt*combeffpst;
                end
               
            end
    
            
       end