%rewritten extensively 4.16.09,
%to incorporate new figure requirements.

%takes sumdyn,
%and a ps (plotstruct).
%ps.minx - minimum x value to include
%ps.maxx - maximum x value to include.
%ps.col,  col={ 'k' 'r','c'}
%ps.colvec
%ps.type, (pct or off or dis)
%ps.addx
%ps.excludebs
%ps.plotavg=1
%ps.comb
%ps.normeff=0

function [outvlaczcomb,outvlmuzcomb]=plotinitdynamics(sumdyn,sumbs,ps)
sdynln=length(sumdyn)
tm_matrix=0:1:ps.maxx
bsindcomb=[];
sumdyncomb=[];
 outstruct=[];   
% 
% % first find the number of traces
% for ii=1:length(sumdyn)
%     if((max(sumdyn(ii).tms))>=ps.minx&ii~=ps.excludeind)
%         bsindcomb=[bsindcomb sumdyn(ii).bsnum];
%         sumdyncomb=[sumdyncomb ii];
%         cnt=cnt+1;
%     end
% end
outvlaczcomb=zeros(sdynln,ps.maxx+1);
outvlmuzcomb=zeros(sdynln,ps.maxx+1);
ct=1;
for ii=1:length(sumdyn)
    
    smcr=sumdyn(ii);
    smcr.tms=[0;smcr.adjtms];
    crsumbs=sumbs(smcr.bsnum);
    crbasruns=crsumbs.basruns(end);
    if(ps.insert)
        acinsert=crsumbs.acz(1,crbasruns);
        muinsert=crsumbs.muz(1,crbasruns);
    else
        acinsert=0
        muinsert=0
    end
   
    aczvl=[acinsert smcr.acz]
    muzvl=[muinsert smcr.muz]
    if(smcr.acz(1)>0)
        smcr.drxn=1;
    else
        smcr.drxn=0;
    end
%     if(exist('ps.addx'))
%         smcr.tms=smcr.tms+addx;
%               
    if(ps.flip)
        if(smcr.drxn==0)
            aczvl=-aczvl
            muzvl=-muzvl
        end
    end
            
    [acvls]=interp1(smcr.tms,aczvl,tm_matrix);
    [muvls]=interp1(smcr.tms,muzvl,tm_matrix);
    
%                 slope=calcslope(smcr.tms(indx),yvl(indx));
%                 equalflag=1;
  
  ln=length(tm_matrix);
    if(max(smcr.tms)>ps.minx)
        outvlaczcomb(ct,1:ln)=acvls(1:ln)
        outvlmuzcomb(ct,1:ln)=muvls(1:ln)
        ct=ct+1;
    
        if(ps.plotraw)

            plot(tm_matrix,outvlaczcomb(ct-1,:),'Color',ps.ac_col,'Linewidth',2)
            hold on;
            
            plot(tm_matrix,acvls,'o','MarkerSize',5,'MarkerFaceColor',ps.ac_col)
            plot(tm_matrix,muvls,'o','MarkerSize',5,'MarkerFaceColor',ps.mu_col)
%             hold on;
            plot(tm_matrix,outvlmuzcomb(ct-1,:),'Color',ps.mu_col,'Linewidth',2)
        end
    end
end
            
 inds=[1:length(outvlaczcomb)]

 [outacmn,outacstd]=calcmeanstder2(outvlaczcomb(inds,:));
[outmumn,outmustd]=calcmeanstder2(outvlmuzcomb(inds,:));

rvacmn=[outacmn]
rvmumn=[outmumn]
rvacer=[outacstd]
rvmuer=[outmustd]


xvec=[0:ps.maxx]
fillx=[xvec xvec(end:-1:1)]

muvec=[rvmumn+rvmuer]
muvec2=[rvmumn-rvmuer]
filly=[muvec muvec2(end:-1:1)]

acvec=[rvacmn+rvacer]
acvec2=[rvacmn-rvacer]
filly2=[acvec acvec2(end:-1:1)]

mufillcol=[0.82 0.82 0.82]
acfillcol=[0.92 0.96 0.98]



% %     mupts=mumean(1:length(indvl))
%     acpts=acmean(1:length(indvl))
%     yvec=[avls.initmean{notevl}*ones(length(mupts),1);mupts(end:-1:1)']
%     

    hold on;
%     fill(xvec,yvec,acfillcol);
% filly([1 end])=0
% filly2([1 end])=0
    fill(fillx,filly,acfillcol,'edgecolor','w');
    fill(fillx,filly2,mufillcol,'edgecolor','w');
plot([xvec], rvacmn,'Color',ps.ac_col,'Linewidth',3)
hold on;
plot([xvec],rvmumn,'Color',ps.mu_col,'Linewidth',3)

if(ps.plotarrow)
    plotarrows2(xvec,outacmn,outmumn,ps);
end
     

% 
% subplot(121)
% plot([1:ps.maxx],outpctmn,'Color','g','Linewidth',4);
% plot([1:ps.maxx;1:ps.maxx],[outpctmn-outpctstd;outpctmn+outpctstd],'g')
% subplot(122)
% plot([1:ps.maxx],outeffmn,'Color','g','Linewidth',4);
% plot([1:ps.maxx;1:ps.maxx],[outeffmn-outeffstd;outeffmn+outeffstd],'g')
% 

% outstruct.effvls=outvleffcomb;
% outstruct.yvls=outvlcomb;
% outstruct.effmn=outeffmn;
% outstruct.effsd=outeffstd;
% outstruct.yvlmn=outpctmn;
% outstruct.yvlsd=outpctstd;
% outstruct.bsvls=bsindcomb;
% outstruct.dynind=sumdyncomb;
% outstruct.yslope=yslope;
% outstruct.effslope=effslope;

function [normout]=normvls(vls,norm)
    if(norm)
        normout=vls./vls(1);
    else
        normout=vls;
    end

function [yvl,yvl2,equalflag]=geteffvl(smcr,indx,ps);
     if(ps.comb&ps.normeff)
          %first normalize each
          normtargeff=norm_onevls(smcr.targeff(indx));
          normctrleff=norm_onevls(smcr.contreff(indx));
          yvl=mean([normtargeff;normctrleff]);
          yvl2=yvl;
          equalflag=1;
     elseif(~ps.comb&~ps.normeff)
            yvl=smcr.targeff(indx);
            yvl2=smcr.contreff(indx);
            equalflag=0;
     elseif(ps.comb&~ps.normeff)
            yvl=mean([smcr.targeff(indx); smcr.contreff(indx)]);
            yvl2=yvl;
            equalflag=1;
     elseif(~ps.comb&ps.normeff)
             yvl=norm_onevls(smcr.targeff(indx));
             yvl2=norm_onevls(smcr.contreff(indx));
             equalflag=0;
     end
                
function [yout]=norm_onevls(yin);
    yout=yin/mean(yin);
    
    
    function [slope]= calcslope(xin,yin);
        
        s=polyfit(xin,yin',1);
        slope=s(1);