%written 11.11.08, meant to compare offset runs across birds, input is a
%ps, a plot struct, and shfts, which is outputted by shiftanal2.
function [ax]=shiftsumplot3(shs)
   
 
 %find which values are targeted and which values are controls.
 for ii=1:length(shs)
  sh=shs(ii);
  rtype=sh.ntype;
  if rtype=='targ'
    numrns=length(sh.subruns)
    clear ax1 ax2
    for jj=1:numrns
        ax1(jj)=subplot(211)
        ax2(jj)=subplot(212);
        
    
    end
          plotpct2(sh, ax1);
        ploteff2(sh,ax2)             
  end
 end