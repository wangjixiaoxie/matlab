%written 9/7/09 to create sumbs struct for stimdata

function [shsall, shsrev, sumbs] = metstimanal(bs, brdind)
global revstdlimit 
global revnumrunslimit
shsall=[];


revstdlimit=1.5;
revnumrunslimit=4;


minstd=0
mineffn=10
%first calculate the shs values for this run
for ct=1:length(brdind)
    bvl=brdind(ct)
    [avls]=loadsumdata2(bs,bvl);
    %shs is a structure for each run
    [shsout,shsrev]=calcshsvals(avls,bs(bvl),minstd,mineffn);
    
    sumbs(bvl).acz=avls.catchz
%     sumbs(bvl).acprez=avls.acprez
%     sumbs(bvl).acpstz=avls.acpstz
%     sumbs(bvl).rawtimes=avls.rawtimes;
%     sumbs(bvl).mulist=avls.mulist;
    sumbs(bvl).initmean=avls.mnbas;
    sumbs(bvl).initsd=avls.stdbas;
    sumbs(bvl).bname=bs(bvl).bname;
    sumbs(bvl).bsnum=bvl;
%     sumbs(bvl).acerrz=avls.acerracz
    sumbs(bvl).acmean=avls.ctmean;
    sumbs(bvl).mumean=avls.fbmean;
    if(isfield(avls,'STAN_RUNS'))
        sumbs(bvl).STANRUNS=avls.STAN_RUNS
    end
%     sumbs(bvl).acn=avls.acn;
    sumbs(bvl).mun=avls.crctind;
%     sumbs(bvl).acstdv=avls.acstdv;
%     sumbs(bvl).mustdv=avls.mustdv;
    sumbs(bvl).muz=avls.stimz
    sumbs(bvl).pct=((avls.catchz-avls.stimz)./avls.catchz)*100;
    sumbs(bvl).initshiftind=bs(bvl).initshiftind;
    sumbs(bvl).basruns=bs(bvl).basruns;
%     sumbs(bvl).muerrz=avls.mustderrz
    
%     for ii=1:length(avls.acstdv(:,1))
%         sumbs(bvl).acstdz(ii,:)=avls.acstdv(ii,:)/avls.initsd{ii}
%         sumbs(bvl).mustdz(ii,:)=avls.mustdv(ii,:)/avls.initsd{ii}
%     end
   
    for ii=1:length(shsout)
        sumbs(bvl).ntind=bs(bvl).ntind;
        sumbs(bvl).sfact=bs(bvl).sfact;
        sumbs(bvl).allnote=shsout(ii).allnote;
%         sumbs(bvl).effvls=shsout(ii).effvls;
%         sumbs(bvl).combeff=shsout(ii).combeff;
        sumbs(bvl).tmvec=shsout(ii).tmvec;
        crtmvec=shsout(ii).tmvec

        sumbs(bvl).basruns=bs(bvl).basruns
        
        [out,sortind]=sort(crtmvec(shsout(ii).shiftruns));
        sumbs(bvl).shiftruns{ii}=shsout(ii).shiftruns(sortind);
        sumbs(bvl).subruns{ii}=shsout(ii).subruns;
        sumbs(bvl).initind{ii}=shsout(ii).revind;
        sumbs(bvl).flr_wnon{ii} =shsout(ii).flr_wnon;
        sumbs(bvl).flr_wnoff{ii}=shsout(ii).flr_wnoff;
        sumbs(bvl).flrtmvec=shsout(ii).flrtmvec;
        
%         sumbs(bvl).pct=shsout(ii).pct;
%         sumbs(bvl).offz=shsout(ii).offz;
        sumbs(bvl).drxn{ii}=shsout(ii).drxn;
        sumbs(bvl).asympvl{ii}=shsout(ii).asympvl;
        sumbs(bvl).asympruns{ii}=shsout(ii).asympind
        sumbs(bvl).muasympdist{ii}=shsout(ii).muasympdist;
        sumbs(bvl).acasympdist{ii}=shsout(ii).acasympdist;
    if(~isempty(shsrev))
        for ii=1:length(shsrev)
         
      
          [out,sortind]=sort(crtmvec(shsrev(ii).runs));
            tmprevruns=shsrev(ii).runs(sortind);
            stanind=ismember(tmprevruns,sumbs(bvl).STANRUNS)
            sumbs(bvl).revruns{ii}=tmprevruns(stanind);
        sumbs(bvl).revflron{ii}=shsrev(ii).flr_wnon;
        sumbs(bvl).drxnrev{ii}=shsrev(ii).drxn;
        end
    else
        sumbs(bvl).revruns=[];
    end
    
    shsall=[shsall shsout];
end

end  
    
%
% %     [shs]=combeffvls(shsall);
%     [sumshs.meaneff,sumshs.stderreff]=calcavevl(shs,'eff',numdays);
%     [sumshs.pct, sumshs.stderrpct]=calcavevl(shs,'pct',numdays);
    sumbs=calcwnoff(sumbs);

    test=2
 

function [shs, shsrev]=calcshsvals(avls,bs,minstd,mineffn)
     wn=avls.wn
     
     if isfield(avls,'wnrev')
        wnrev=avls.wnrev;
     else
         wnrev=[];
         shsrev=[]
     end
     ntind=bs.ntind;
     contrind=bs.contrind;
     allind1=[ntind contrind]     
     contrvl=[zeros(length(ntind),1); ones(length(contrind),1)]
     
     clear subind runs diff
     outruns=[];
%      for zz=1 length(allind1)
    for zz=1
         ii=allind1(zz);
        %each row (jj) is a different shift
        %this loop sets the runs{1}(jj,:)-lim runs, 
        %and runs{2}(jj,:)-all runs
        %also sets the mxshift(jj), drxn{jj}
        for jj=1:length(wn(ii).on(:,1))
                
            [flr_wnon, flr_wnoff, subind,flrtmvec,tmvec]=calcshiftvls(wn, ii, jj,avls)
            [outruns,maxshift(jj),drxn{jj},asympvl{jj},muasympdist,acasympdist]=findruns(avls,subind,ii,jj,zz,flrtmvec,contrvl,minstd,bs);
            [shs(jj)]=calcoutvls(flrtmvec, flr_wnon, flr_wnoff, avls, maxshift(jj), outruns, allind1,ii,mineffn,drxn{jj},asympvl{jj},muasympdist,acasympdist,tmvec);
        end

        if(~isempty(wnrev))
            if(~isempty(wnrev(ii).on))
                for kk=1:length(wnrev(ii).on(:,1))
                    [flr_wnon,flr_wnonff,subind,flrtmvec]=calcrevshiftvls(wnrev,ii,kk,avls);
                    [revruns,drxn{kk}]=findrevruns(avls,subind,ii,kk,zz,flrtmvec);
                    [shsrev(kk)]=calcoutrevvls(flrtmvec,flr_wnon,flr_wnoff,avls,revruns,allind1,ii,drxn{kk});
                end
            else
                shsrev=[];
            end
        else
            shsrev=[];
        end
    end

    function [outruns, mxshift, drxn,asympvl,muasympdist,acasympdist]=findruns(avls,subind,ii,jj,zz,flrtmvec,contrvl,minstd,bs);
        avls.acz=avls.catchz
        avls.muz=avls.stimz
        
        for subindvl=1:3
                    runs=subind{subindvl};
%                     matchind=find(ismember(crind,avls.mulist));
%                     runs=crind(matchind);
                    [sortruns,sortind]=sort(flrtmvec(runs));
                    %only reset outruns if run is a target run
                    if(sortind&contrvl(zz)==0)
                        outruns{subindvl}=runs(sortind);
                        outrnstmp=runs(sortind);
                    
                        if(subindvl==1)
                            aczvls=avls.acz(ii,outrnstmp);
                            absvls=abs(aczvls);
                            indabs=find(absvls>minstd);
                            outruns{subindvl}=outrnstmp(indabs);
                        end
                    end
        end
                aczall=avls.acz(ii,outruns{2});
                [tmp,mxind]=max(abs(aczall));
                mxshift=aczall(mxind)     
                if(mxshift>0)
                    drxn='up'
                else
                    drxn='do'
                end 
                
                 %calculate outruns{4}, asympind.
                
                    [outruns{4},asympvl,muasympdist,acasympdist]=calcasympind(avls.acz(ii,:), avls.muz(ii,:),outruns{2}, mxshift,bs)
                
                test=1
        
        function [revruns,drxn]=findrevruns(avls,subind,ii,jj,zz,flrtmvec);
%                     
                avls.acz=avls.catchz
                avls.muz=avls.stimz
%                 matchind=find(ismember(subind,avls.mulist));
                    runs=subind;
                    
                    
                    [sortruns,sortind]=sort(flrtmvec(runs));
                    %only reset outruns if run is a target run
                        outrnstmp=runs(sortind);
                        if(~isempty(outrnstmp))
                            if(avls.acz(ii,outrnstmp(1))>0)
                                drxn='up'
                            else
                                drxn='do'
                            end
                        else
                            drxn='up'
                        end
                        revruns=outrnstmp;
                
               
                    
                
        function [flr_wnon, flr_wnoff, subind,flrtmvec,tmvec]=calcshiftvls(wn,ntin,shin,avls)

                flr_wnonall(shin,:)=floor(wn(ntin).on(shin,:))
                flr_wnoffall(shin,:)=floor(wn(ntin).off(shin,:));
                indon=find(flr_wnonall(shin,:));
                indoff=find(flr_wnoffall(shin,:))
                flr_wnon=flr_wnonall(shin,indon);
                flr_wnoff=flr_wnoffall(shin,indoff);
                tmvec=avls.adjtimes(:,1);
                flrtmvec=floor(avls.adjtimes(:,1));
            %subind{1} are all the inds, subind{2} are the asymp inds.
                subind{1}=find(flrtmvec>=flr_wnon(end)&flrtmvec<flr_wnoff(end));    
                subind{2}=find(flrtmvec>=flr_wnon(1)&flrtmvec<flr_wnoff(end));    
                subind{3}=1:length(flrtmvec);
                
          function [flr_wnon, flr_wnoff, subind,flrtmvec]=calcrevshiftvls(wnrev,ntin,shin,avls)

                flr_wnonall(shin,:)=floor(wnrev(ntin).on(shin,:))
                flr_wnoffall(shin,:)=floor(wnrev(ntin).off(shin,:));
                indon=find(flr_wnonall(shin,:));
                indoff=find(flr_wnoffall(shin,:))
                flr_wnon=flr_wnonall(shin,indon);
                flr_wnoff=flr_wnoffall(shin,indoff);
                
                flrtmvec=floor(avls.adjtimes(:,1));
            %subind{1} are all the inds, subind{2} are the asymp inds.
                subind=find(flrtmvec>=flr_wnon(end)&flrtmvec<flr_wnoff(end));        

function  [shs]=calcoutvls(flrtmvec, flron, flroff, avls,mxshift,outruns,allind1,ntind,mineffn,drxn,asympvl,muasympdist,acasympdist,tmvec)               
                avls.acz=avls.catchz
                avls.muz=avls.stimz
                shs.flrtmvec=flrtmvec;
                shs.tmvec=tmvec;
                shs.flr_wnon=flron;
                shs.drxn=drxn;
                shs.asympvl=asympvl;
  %to calculate this/need to switch into mu coordinates.
                shs.flr_wnoff=flroff;
                shs.acz=avls.acz(ntind,:);
                shs.muz=avls.muz(ntind,:);
%                 shs.acerrz=avls.acerracz(ntind,:)
                shs.mxshift=mxshift
%                 shs.muerrz=avls.mustderrz(ntind,:)
                shs.shiftruns=outruns{2}
                shs.subruns=outruns{1};
                shs.allruns=outruns{3}(1,:);
                shs.asympind=outruns{4};
                shs.muasympdist=muasympdist;
                shs.acasympdist=acasympdist;
                shs.allnote=allind1;
%                 shs.bname=avls.bname'
                %if control note
                shs.ntype='targ'
                 %this is for all vals
%                  [shs.effvls,shs.combeff,shs.ac_cv, shs.mu_cv]=get_effvls(avls,allind1,mineffn);
%                  [shs.offz,shs.offerrz,shs.off,shs.offerr]=getoff(avls,allind1);
                 %this will selected indices for reversion runs
                 %and selected off values, pct values, and errors.
                 [shs.revind]=select_rev(shs,avls);
                 numcontrind=length(allind1)-1;
                %currently just for subruns vals, but I am going to change this. 
%                 [shs.pct,shs.pcter]=getpctvls2(avls,allind1(1),shs);
%                 
                
  
  function  [shsrv]=calcoutrevvls(tmvec, flron, flroff, avls,outruns,allind1,ntind,drxn)               
                shsrv.drxn=drxn;
                shsrv.flrtmvec=tmvec;
                shsrv.flr_wnon=flron;
  %to calculate this/need to switch into mu coordinates.
                shsrv.flr_wnoff=flroff;
                shsrv.runs=outruns;
                
%                 shsrv.bname=avls.bname'
                %if control note
                %currently just for subruns vals, but I am going to change this. 
                   

                
                
  %the point of the function is to go throuch each of the shs subruns, combine the multiple
  %effvls into a daily estimate, and then interpolate the estimate
%   and interpolate through all the gap days                      
  function [mnvl,stdvl]=calcavevl(shs,yfieldname,numdays)
      initval=-1000
      %initialize effinterpvl matrix to 0
      %number of rows is number of shifts
      %ncol is numdays.
          
  %FIX THIS HERE...CREATING TOO MANY EFFINTERPVLS OUTS!!
  %create one output for run
  %and then combine them
  interpcomb=[];
  for ii=1:length(shs)
          shscr=shs(ii);
          if ~isempty(shs(ii).subruns)
                [dys]=calcxvls(shscr)
          else
              dys=[];
          end 
              nrow=length(shs(ii).subruns)
                ncol=numdays;
          
          interpvlstmp=zeros(ncol,1);
          interpvlstmp=interpvlstmp-1000;
          if ~isempty(dys)
            
            if(yfieldname=='eff')
                cmd=['yvls=shscr.combeff(shscr.subruns)'];
                eval(cmd);
            else
                cmd=['yvls=shscr.pct'];
                eval(cmd);
            end
             
          [interpvlstmp(dys(1):1:dys(end))]=interpvls(yvls,dys);
          end
            interpcomb=[interpcomb;interpvlstmp(1:numdays)']  
      end
          %need to calculate the mean column by column excluding zero
         
      
           [mnvl,stdvl]=calcmatmean(interpcomb,numdays,initval);  
            %loop through days and calculate the difference with previous
            %day
            %if diff>1, interpola
    %returns shs.combeff 
   function [shs]=combeffvls(shs);              
        for shin=1:length(shs)
            shscr=shs(shin);
            allnote=shscr.allnote;
            nmnote=length(allnote);
            shs(shin).combeff=mean(shscr.effvls(allnote,:),1)
            if (nmnote>1)
                shs(shin).stderref=std(shscr.effvls(allnote,:),1)/sqrt(nmnote);
            else
                shs(shin).stderref=0;
            end
        end
                
   function [dys]=calcxvls(shscr) 
                subrns=shscr.subruns
                flrtmvc=shscr.flrtmvec;
                ofst_tm=flrtmvc-flrtmvc(subrns(1));
                %this is to avoid any zeros
                dys=ofst_tm(subrns)+1;
   
   

       %loop through days and calculate the difference with previous
            %day
            %if diff>1, interpolate  
            %DO NOT INCLUDE -1000 values in mean
    function [mnvl,stderrvl]=calcmatmean(initmat,numdays,initval);  
        for ii=1:numdays
           initmatcol=initmat(:,ii);
           ind=find(initmatcol~=initval);
           mnvl(ii)=mean(initmatcol(ind));
     stderrvl(ii)=std(initmatcol(ind))/sqrt(length(ind));
        end

    %write this function
    %
    function [ydiffz,ydifferrz,ydiff,ydifferr]=getoff(avls,allind1);    
        for ii=1:length(allind1)
            ntind=allind1(ii);
            acz=avls.acz(ntind,:);
            muz=avls.muz(ntind,:)
            acerrz=avls.acerracz(ntind,:);
            muerrz=avls.mustderrz(ntind,:);
            
            acmean=avls.acmean(ntind,:);
            mumean=avls.mumean(ntind,:)
            acerr=avls.acstderr(ntind,:);
            mustderr=avls.mustderr(ntind,:);
        
            ydiffz(ntind,:)=-acz+muz;
            ydifferrz(ntind,:)=(acerrz+muerrz)/2;
            
            ydiff(ntind,:)=mumean-acmean;
            ydifferr(ntind,:)=(acerr+mustderr)/2;
            
        end
        
        
    %this will selected indices for reversion runs
    %and selected off values, pct values, and errors.
    function[shsrevind]=select_rev(shs,avls);
       global revstdlimit
       global revnumrunslimit
       
        shiftind=shs.shiftruns;
        acz=avls.acz
        ntind=shs.allnote(1);
        absacz=abs(acz(ntind,shiftind));
        ind=find(absacz>revstdlimit);
        if(length(ind)>1)
            if length(ind)>revnumrunslimit
                ind=ind(1:revnumrunslimit);
            else
                ind=ind;
            end
                shsrevind=shiftind(ind);
          
        elseif ind==1
            shsrevind=shiftind(ind(1));
                
        else
            shsrevind=[]
        end

   function [sumbs]=calcwnoff(sumbs)
        
       for ii=1:length(sumbs)
         shiftrns=sumbs(ii).shiftruns;
         
         revruns=sumbs(ii).revruns;
         asympruns=sumbs(ii).asympruns;
         ntind=sumbs(ii).ntind;   
         for jj=1:length(shiftrns)
                runs=shiftrns{jj};
                tms=sumbs(ii).flrtmvec(runs);
                [out_tms]=calc_exact_tms(sumbs(ii).tmvec(runs));
               
                
                
                
                baswntime=sumbs(ii).flr_wnon{jj}(1);
                %the plus one is so that if inactivation is same day as wn,
                %then that is listed as day 1.
                sumbs(ii).adjshifttms{jj}=tms-baswntime;
                sumbs(ii).exactshifttms{jj}=out_tms-baswntime;
                sumbs(ii).baswntime{jj}=baswntime;
               
         end
         for jj=1:length(revruns)
             runs=revruns{jj};
             tms=sumbs(ii).flrtmvec(runs);
             [out_tms]=calc_exact_tms(sumbs(ii).tmvec(runs));
             baswntime=sumbs(ii).revflron{jj};
             sumbs(ii).revshifttms{jj}=tms-baswntime+1;
              sumbs(ii).exactrevtms{jj}=out_tms-baswntime+1;
         end
        if(isfield(sumbs(ii),'asympruns'))
            if(~isempty(sumbs(ii).asympruns))
                for jj=1:length(asympruns)
            runs=asympruns{jj}
            tms=sumbs(ii).flrtmvec(runs);
            [out_tms]=calc_exact_tms(sumbs(ii).tmvec(runs));
                baswntime=sumbs(ii).flr_wnon{jj}(1);
                %the plus one is so that if inactivation is same day as wn,
                %then that is listed as day 1.
                sumbs(ii).asympshifttms{jj}=tms-baswntime;
                sumbs(ii).exactasymptms{jj}=out_tms-baswntime;
                end
            end
         end 


       end
  
       
 function [recenter_tms]=calc_exact_tms(intms)
      
                modtms=mod(intms,floor(intms));
                recenter_tms=(modtms*24-7)/14+floor(intms);
     

     function [sumbs]=overalleff(sumbs)  
       for ii=1:length(sumbs)
         crbs=sumbs(ii);
         
             ind=find(crbs.combeff>0&isnan(crbs.combeff)==0)
             sumbs(ii).mneff=mean(crbs.combeff(ind));
             sumbs(ii).stdeff=std(crbs.combeff(ind))
             
        end

     function [asympind,asympvl,muasympdist,acasympdist]=calcasympind(aczvls, muzvls,ind, mxshift,bs)
        
         
         %sdnet is how many standard deviations out to calculate the mean.
         sdnet=1.5;
         %sd pass is distance from mean to include as part of asympvks
         sdpass=1
         
         
            if (mxshift<0)
                vls=-aczvls(ind)
                mxshift=-mxshift;
            else
                vls=aczvls(ind)
            end
            
            %HACK
%              if(bs.skipasymp)
%                  vls=aczvls
%                  ind=1:length(aczvls)
%                  mxshift=vls
%              end
         
         
          if isfield(bs,'asympvl')
             mnvl=bs.asympvl
             indcheck=1:length(aczvls(ind))
          else
            indcheck=find(abs(mxshift-vls)<sdnet)
            mnvl=mean(vls(indcheck))
          end
         
         indpass=find(abs(vls(indcheck)-mnvl)<sdpass&abs(vls(indcheck))>1.5);
        
         asympind=ind(indcheck(indpass));
         asympvl=mean(aczvls(ind(indcheck(indpass))))
         muasympdist=asympvl-muzvls(ind(indcheck(indpass)));
         acasympdist=asympvl-aczvls(ind(indcheck(indpass)));
         if(asympvl<0)
             muasympdist=-muasympdist
             acasympdist=-acasympdist
         end
         
         