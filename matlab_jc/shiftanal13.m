%5.4.10
%modified to calculate cnscasympind, call outside function calcasympind




%currently written to pick out asymptotic runs
% %with >2sd shifts, and same wn contingencies.
% %rewritten 11.7.08 in order to calculate for all muruns these values
% %and output shfts(ntcnt).acz
%             shfts.muz
%             shfts.offset
%             shfts.acerrz
%             shfts.muerrz
%             shfts.diffz
%             shfts.differrz
%             shfts.bs
%             shfts.nshfts
%             shfts.shftind -- this is a series of inds for each shift
%             shfts.targ=1  0 is control.
%             shfts.subind
%             shfts.basind
%             shfts.ntind
%[shft_s]=shiftanal2(bs,brdind)

%minstd threshold for subruns

function [sumbs,shsall,sumshs, shsrev] = shiftanal12(bs, brdind)
global revstdlimit 
global revnumrunslimit

revstdlimit=1.5;
revnumrunslimit=4;

shsall=[];
minstd=0
outct=1;
numdays=6;
mineffn=10;
MINCVLEN=5;
%first calculate the shs values for this run
for ct=1:length(brdind)
    bvl=brdind(ct)
    [avls,graphvals]=loadsumdata2(bs,bvl);
    [allmean,allz]=calcallz(avls);
   
    avls.allz=allz
    %shs is a structure for each run
    [shsout,shsrev]=calcshsvals(avls,graphvals,bs(bvl),minstd,mineffn,MINCVLEN);
    
    sumbs(bvl).acz=avls.acz
%     sumbs(bvl).acprez=avls.acprez
%     sumbs(bvl).acpstz=avls.acpstz
    sumbs(bvl).tmvc=avls.rawtimes;
    sumbs(bvl).mulist=avls.mulist;
    sumbs(bvl).aclist=avls.aclist;
    
    sumbs(bvl).bname=avls.bname;
    sumbs(bvl).acerrz=avls.acerracz
    sumbs(bvl).acmean=avls.acmean;
    sumbs(bvl).mumean=avls.mumean;
    sumbs(bvl).allmean=allmean
    sumbs(bvl).acn=avls.acn;
    sumbs(bvl).mun=avls.mun;
    sumbs(bvl).acstdv=avls.acstdv;
    sumbs(bvl).mustdv=avls.mustdv;
    sumbs(bvl).muz=avls.muz
    
    sumbs(bvl).muerrz=avls.mustderrz
    
    for ii=1:length(avls.acstdv(:,1))
        sumbs(bvl).acstdz(ii,:)=avls.acstdv(ii,:)/avls.initsd{ii}
        sumbs(bvl).mustdz(ii,:)=avls.mustdv(ii,:)/avls.initsd{ii}
    end
    
   
    for ii=1:length(shsout)
        sumbs(bvl).ntind=bs(bvl).ntind;
        sumbs(bvl).sfact=bs(bvl).sfact;
        sumbs(bvl).allnote=shsout(ii).allnote;
        sumbs(bvl).effvls=shsout(ii).effvls;
        sumbs(bvl).combeff=shsout(ii).combeff;
        sumbs(bvl).allcv=shsout(ii).allcv;
        sumbs(bvl).basruns=bs(bvl).basruns
        sumbs(bvl).shiftruns{ii}=shsout(ii).shiftruns;
        sumbs(bvl).allshiftruns{ii}=shsout(ii).allshiftruns;
        sumbs(bvl).limshiftruns{ii}=shsout(ii).limshiftruns;
         sumbs(bvl).ac_cv=shsout(ii).cv;
        sumbs(bvl).mu_cv=shsout(ii).cv;
        sumbs(bvl).subruns{ii}=shsout(ii).subruns;
        sumbs(bvl).initind{ii}=shsout(ii).revind;
        sumbs(bvl).flr_wnon{ii} =shsout(ii).flr_wnon;
        sumbs(bvl).flr_wnoff{ii}=shsout(ii).flr_wnoff;
        sumbs(bvl).flrtmvec=shsout(ii).flrtmvec;
        sumbs(bvl).pct=shsout(ii).pct;
        sumbs(bvl).offz=shsout(ii).offz;
        sumbs(bvl).drxn{ii}=shsout(ii).drxn;
        sumbs(bvl).asympvl{ii}=shsout(ii).asympvl;
        sumbs(bvl).asympruns{ii}=shsout(ii).asympind
        sumbs(bvl).cnscasympruns{ii}=shsout(ii).cnscasympruns;
        for jj=1:length(sumbs(bvl).allnote)
            crnote=sumbs(bvl).allnote(jj);
            sumbs(bvl).initmean(crnote)=avls.initmean{crnote};
            sumbs(bvl).initsd(crnote)=avls.initsd{crnote};
        end
%         sumbs(bvl).muasympdist{ii}=shsout(ii).muasympdist;
%         sumbs(bvl).acasympdist{ii}=shsout(ii).acasympdist;
    shsall=[shsall shsout];
    end
for ii=1:length(shsrev)
        sumbs(bvl).revruns{ii}=shsrev(ii).runs
        sumbs(bvl).revasympruns{ii}=shsrev(ii).asympruns;
        sumbs(bvl).cnscrevasympruns{ii}=shsrev(ii).cnscasympruns;
        sumbs(bvl).revasympvl{ii}=shsrev(ii).asympvl;
        sumbs(bvl).allrevruns{ii}=shsrev(ii).allruns
        sumbs(bvl).revflron{ii}=shsrev(ii).flr_wnon;
        sumbs(bvl).drxnrev{ii}=shsrev(ii).drxn;
    end
    
    
   
  [shiftlist,asympshiftlist,revlist]=mklists(sumbs(bvl));
  sumbs(bvl).shiftlist=shiftlist;
  sumbs(bvl).asympshiftlist=asympshiftlist;
  sumbs(bvl).revlist=revlist;
  sepbasruns=mksepbasruns(sumbs(bvl));
  sumbs(bvl).sepbasruns=sepbasruns;
  if(bvl==6)
      sumbs(bvl).sepbasruns{1}=24;
  end
  
end  
    
%
    [shs]=combeffvls(shsall);
    [sumshs.meaneff,sumshs.stderreff]=calcavevl(shs,'eff',numdays);
    [sumshs.pct, sumshs.stderrpct]=calcavevl(shs,'pct',numdays);
    sumbs=calcwnoff(sumbs);
    
    test=2

    function[sepbasruns]=mksepbasruns(crbs)
        for ii=1:length(crbs.allshiftruns)
            
            [mintm,maxtm]=get_edgetms(crbs,ii);
            bastms=crbs.flrtmvec(crbs.basruns);
            ind=find(bastms>mintm&bastms<maxtm)
            if(~isempty(ind))
                sepbasruns{ii}=crbs.basruns(ind);
            else
                sepbasruns{ii}=[];
            end
            
        end
function[crmean,crz]=calcallz(avls)
   
    for ii=1:length(avls.adjvls)
        for jj=1:length(avls.adjvls{ii})
            crmean(jj,ii)=mean(avls.adjvls{ii}{jj}(:,2));
            crz(jj,ii)=(crmean(jj,ii)-avls.initmean{ii})./avls.initsd{ii};
        end
    end
    
   

function[mintm,maxtm]=get_edgetms(crbs,indnum)
    cr_runs=crbs.allshiftruns{indnum};
    cr_tms=crbs.flrtmvec(cr_runs);
    maxtm=min(cr_tms);
    
    if(indnum==1)
        mintm=0;
    else
        mintm=crbs.flrtmvec(crbs.allshiftruns{indnum-1}(end));
    end
        
    
function [shlist,asympshlist,revlist]=mklists(crbs)
    shlist=[]
    asympshlist=[];
    revlist=[];
    asymprevlist=[];
    baslist=[];
    
    for ii=1:length(crbs.allshiftruns)
        shlist=[shlist makerow(crbs.limshiftruns{ii})]; 
    end
%     for ii=1:length(crbs.allasympruns)
%         asympshlist=[asympshlist makerow(crbs.allasympruns{ii})];
%     end
if(isfield(crbs,'allrevruns'))
    for ii=1:length(crbs.allrevruns)
        revlist=[revlist makerow(crbs.allrevruns{ii})];
        
    end
end
    
            
        
    
function [shs, shsrev]=calcshsvals(avls,graphvals,bs,minstd,mineffn,MINCVLEN)
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
                
            [flr_wnon, flr_wnoff, subind,flrtmvec]=calcshiftvls(wn, ii, jj,avls)
            [outruns,outrunsall,maxshift(jj),drxn{jj},asympvl{jj}]=findruns(avls,subind,ii,jj,zz,flrtmvec,contrvl,minstd,bs,flr_wnon,flr_wnoff);
            [shs(jj)]=calcoutvls(flrtmvec, flr_wnon, flr_wnoff, avls, maxshift(jj), outruns, outrunsall,allind1,ii,mineffn,drxn{jj},asympvl{jj},MINCVLEN);
        end

        if(~isempty(wnrev))
            if(~isempty(wnrev(ii).on))
                for kk=1:length(wnrev(ii).on(:,1))
                    [flr_wnon,flr_wnonff,subind,flrtmvec]=calcrevshiftvls(wnrev,ii,kk,avls);
                    [revruns,revrunsall,drxn{kk},asympruns{kk},cnscasympruns{kk},asympvl{kk}]=findrevruns(avls,subind,ii,flrtmvec,bs,flr_wnon,flr_wnoff);
                    [shsrev(kk)]=calcoutrevvls(flrtmvec,flr_wnon,flr_wnoff,avls,revruns,revrunsall,allind1,ii,drxn{kk},asympruns{kk},cnscasympruns{kk},asympvl{kk});
                end
            else
                shsrev=[];
            end
        else
            shsrev=[];
        end
    end

    
    %subind{1} ??, subind{2} are the shift
            %inds. subind{3} are all the inds
    function [outruns, outrunsall,mxshift, drxn,asympvl]=findruns(avls,subind,ii,jj,zz,flrtmvec,contrvl,minstd,bs,flrwnon,flrwnoff);
        for subindvl=1:3
                    crind=subind{subindvl};
                    
                    
                    matchind=find(ismember(crind,avls.mulist));
                    runs=crind(matchind);
                    allruns=crind;
                    [sortruns,sortind]=sort(flrtmvec(runs));
                    [sortrunsall,sortindall]=sort(flrtmvec(allruns));
                    
                    %only reset outruns if run is a target run
                    if(sortind&contrvl(zz)==0)
                        outruns{subindvl}=runs(sortind);
                        outrunsall{subindvl}=allruns(sortindall);
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
                mxshift=aczall(mxind);
               
               
                
                if(mxshift>0)
                    drxn='up'
                else
                    drxn='do'
                end 
                
                 %calculate outruns{4}, asympind.
                 ps.ntind=bs.ntind;
                 ps.indin=outrunsall{2};
                 ps.flrwnon=flrwnon;
                 ps.flrwnoff=flrwnoff;
                 ps.flrtmvec=flrtmvec;
                 ps.type='phar'
                 aczvls=avls.acz(ii,:)
                [outruns{4},outruns{5},asympvl]=calcasympind(avls,avls.allz(:,ii),ps)
                
                
                %                 muind=find(ismember(avls.mulist,outruns{4}));
%                 
%                 
%                 preind=avls.aclist(muind,1);
%                 postind=find(avls.aclist(muind,2)
%                 test=1
        
        function [revruns,revrunsall,drxn,asympruns,cnscasympruns,asympvl]=findrevruns(avls,subind,ii,flrtmvec,bs,flrwnon,flrwnoff)
                    matchind=find(ismember(subind,avls.mulist));
                    runs=subind(matchind);
                    runsall=subind;
                    
                    [sortruns,sortind]=sort(flrtmvec(runs));
                    [sortrunsall,sortindall]=sort(flrtmvec(runsall));
                    %only reset outruns if run is a target run
                        outrnstmp=runs(sortind);
                        outrunstmpall=runsall(sortindall);
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
                       if(exist('outrunstmpall'))
                        revrunsall=outrunstmpall;
                       else
                           revrunsall=[];
                       end
               [asympruns,cnscasympruns,asympvl]=calcrevasympind(avls.allz(:,ii), avls.muz(ii,:),revruns,bs,flrwnon,flrwnoff,flrtmvec)
                    
                
        function [flr_wnon, flr_wnoff, subind,flrtmvec]=calcshiftvls(wn,ntin,shin,avls)

                flr_wnonall(shin,:)=floor(wn(ntin).on(shin,:))
                flr_wnoffall(shin,:)=floor(wn(ntin).off(shin,:));
                indon=find(flr_wnonall(shin,:));
                indoff=find(flr_wnoffall(shin,:))
                flr_wnon=flr_wnonall(shin,indon);
                flr_wnoff=flr_wnoffall(shin,indoff);
                
                flrtmvec=floor(avls.adjtimes(:,1));
            %subind{1} are all the asymp inds, subind{2} are the shift
            %inds. subind{3} are all the inds
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

                
                  %subind{4} are all the asymp inds, subind{2} are the shift
            %inds. subind{3} are all the inds
function  [shs]=calcoutvls(tmvec, flron, flroff, avls,mxshift,outruns,outrunsall,allind1,ntind,mineffn,drxn,asympvl,MINCVLEN)               
         
                shs.flrtmvec=tmvec;
                shs.flr_wnon=flron;
                shs.drxn=drxn;
                shs.asympvl=asympvl;
  %to calculate this/need to switch into mu coordinates.
                shs.flr_wnoff=flroff;
                shs.acz=avls.acz(ntind,:);
                shs.muz=avls.muz(ntind,:);
                
                shs.acerrz=avls.acerracz(ntind,:)
                shs.mxshift=mxshift
                shs.muerrz=avls.mustderrz(ntind,:)
                shs.shiftruns=outruns{2}
                %here you need to calculate shiftruns up until the first
                %asymprun; to get all runs find all runs which have a data 
                %less than this max day.
                shs.asympind=outruns{4};
                 shs.allshiftruns=outrunsall{2};
                 for ii=1:length(avls.adjvls)
                    
                    
                        for jj=1:length(avls.adjvls{ii})   
                            if(~isempty(avls.adjvls{ii}{jj}))
                                crptvls=avls.adjvls{ii}{jj}(:,2);
                                shs.cv(ii,jj)=std(crptvls)./mean(crptvls);
                                if(length(crptvls)<MINCVLEN)
                                    
%                                     shs.cv(ii,jj)=NaN;
                                end
                                
%                                 shs.cv(ii,jj)=std(crptvls)./mean(crptvls);
                            end
                        end
                 end
                 
                 
                if(~isempty(outruns{4}))
                    asympinit=shs.asympind(1);
%                     limind=find(shs.shiftruns==asympinit);
%                     limrun=shs.shiftruns(limind);
                    ind=find(avls.tmvc(:,1)<ceil(avls.tmvc(shs.asympind(1),1)));
                    shs.limshiftruns=intersect(shs.allshiftruns,ind);
                else
                    shs.limshiftruns=shs.allshiftruns
                end
                
                
                
                shs.cnscasympruns=outruns{5};
%                 shs.allruns=outrunsall{3};
                shs.subruns=outruns{1};
                shs.allruns=outruns{3}(1,:);
                
%                 shs.muasympdist=muasympdist;
%                 shs.acasympdist=acasympdist;
                shs.allnote=allind1;
                shs.bname=avls.bname'
                %if control note
                shs.ntype='targ'
                 %this is for all vals
                 if(~isempty(allind1))
                    [shs.effvls,shs.combeff,shs.ac_cv, shs.mu_cv,shs.allcv]=get_effvls(avls,allind1,mineffn);
                    [shs.offz,shs.offerrz,shs.off,shs.offerr]=getoff(avls,allind1);
                 %this will selected indices for reversion runs
                 %and selected off values, pct values, and errors.
                    [shs.revind]=select_rev(shs,avls);
                    numcontrind=length(allind1)-1;
                %currently just for subruns vals, but I am going to change this. 
                    [shs.pct,shs.pcter]=getpctvls2(avls,allind1(1),shs);
                 end
                
  
  function  [shsrv]=calcoutrevvls(tmvec, flron, flroff, avls,outruns,outrunsall,allind1,ntind,drxn,asympruns,cnscasympruns,asympvl)               
                shsrv.drxn=drxn;
                shsrv.flrtmvec=tmvec;
                shsrv.flr_wnon=flron;
  %to calculate this/need to switch into mu coordinates.
                shsrv.flr_wnoff=flroff;
                shsrv.runs=outruns;
                shsrv.asympruns=asympruns;
                shsrv.cnscasympruns=cnscasympruns;
                shsrv.asympvl=asympvl;
                shsrv.allruns=outrunsall
                shsrv.bname=avls.bname'
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
            
                ind=1:length(muz);
                
                ydiffz(ntind,:)=-acz(ind)+muz(ind);
                ydifferrz(ntind,:)=(acerrz(ind)+muerrz(ind))/2;
                ydiff(ntind,:)=mumean(ind)-acmean(ind);
                ydifferr(ntind,:)=(acerr(ind)+mustderr(ind))/2;
    
            
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
         if(isfield(sumbs,'revruns'))
            revruns=sumbs(ii).revruns;
         end
         asympruns=sumbs(ii).asympruns;
         ntind=sumbs(ii).ntind;   
         for jj=1:length(shiftrns)
                runs=shiftrns{jj};
                tms=sumbs(ii).flrtmvec(runs);
                baswntime=sumbs(ii).flr_wnon{jj}(1);
                %the plus one is so that if inactivation is same day as wn,
                %then that is listed as day 1.
                sumbs(ii).adjshifttms{jj}=tms-baswntime+1;
               
         end
         if(exist('revruns'))
             if(~isempty(revruns))
                for jj=1:length(revruns)
                    runs=revruns{jj};
                    tms=sumbs(ii).flrtmvec(runs);
                    baswntime=sumbs(ii).revflron{jj};
                    sumbs(ii).revshifttms{jj}=tms-baswntime+1;
                end
             end
         end
        if(isfield(sumbs(ii),'asympruns'))
            if(~isempty(sumbs(ii).asympruns))
                for jj=1:length(asympruns)
            runs=asympruns{jj}
            tms=sumbs(ii).flrtmvec(runs);
                baswntime=sumbs(ii).flr_wnon{jj}(1);
                %the plus one is so that if inactivation is same day as wn,
                %then that is listed as day 1.
                sumbs(ii).asympshifttms{jj}=tms-baswntime+1;
                end
            end
                end 


       end
        
 function [sumbs]=overalleff(sumbs)  
       for ii=1:length(sumbs)
         crbs=sumbs(ii);
         
             ind=find(crbs.combeff>0&isnan(crbs.combeff)==0)
             sumbs(ii).mneff=mean(crbs.combeff(ind));
             sumbs(ii).stdeff=std(crbs.combeff(ind))
             
        end

%      function [asympind,asympindall,asympvl]=calcasympind(avls,indall,ntind)
%         
%          
%          %new algorithm written 2.10.10
%          sd_dist=1.2 ;
%          minshift=1.5;
%          
%          initacind=find(ismember(avls.aclist(:,1),indall));
%          initind=avls.aclist(initacind,1);
%          for ii=1:length(initind)
%              crind=initind(ii);
%              mnvls(ii)=mean(avls.adjvls{ntind}{crind}(:,2));
%          end
%          
%          sdlist=(mnvls-avls.initmean{ntind})/avls.initsd{ntind};
%          
%          
%          sd_difflist=sdlist(2:end)-sdlist(1:end-1);
% %          tmdifflist=tmlist(2:end)-tmlist(1:end-1);
%          
%          
%          ind_belowdiff=find(abs(sd_difflist)<sd_dist);
%          
%          %get the last value for which there is a pair.
%          if(~isempty(ind_belowdiff))
%             lstasympind=ind_belowdiff(end)+1
%             [sortout,ind]=sort(abs(sdlist));
%             asympvl=mean([sdlist(ind(end)) sdlist(ind(end-1))])
% %             [asympvl,ind]=max(abs(sdlist));
% %             asympvl=sdlist(ind);
%             allasympind=find(abs(asympvl-sdlist)<sd_dist&abs(sdlist)>minshift);
% %             mxvl=max(abs(sdlist(allasympind)));
% %             limind=find(abs(mxvl)-abs(sdlist)<sdmax);
% %             allasympind=intersect(allasympind,limind);
%             allasympind=initind(allasympind);
%             init_tm=floor(avls.tmvc(allasympind(1),1));
%             end_tm=ceil(avls.tmvc(allasympind(end),1));
%             
%             %asympindall are all runs within these time bnds.
%             asympindall=find(avls.tmvc(:,1)>=init_tm&avls.tmvc(:,1)<=end_tm);
%             [sortout,ind]=sort(avls.tmvc(asympindall,1));
%             asympindall=asympindall(ind);
%             
%             %so asympind are all mu runs within these time bnds
%             asympind=find(ismember(avls.mulist,asympindall));
%             
%             %asympind are within TKTK
%             
%             asympind=avls.mulist(asympind);
%          
%             
%             
%             %asympindall are all runs within these time bnds
%             
%             
%          else
%              asympind=[];
%              asympindall=[];
%              asympvl=sdlist(end);
%          end
%          
%          
% %          if (mxshift<0)
% %                 vls=-aczvls(ind)
% % %                 vlsall=-aczvls(indall);
% %                 mxshift=-mxshift;
% %             else
% %                 vls=aczvls(ind)
% % %                 vlsall=aczvls(indall);
% %             end
% %          
% %          
% %           if ~isempty(bs.asympvl)
% %              mnvl=bs.asympvl
% %              indcheck=1:length(aczvls(ind))
% %           else
% %             indcheck=find(abs(mxshift-vls)<sdnet)
% % %             indcheckall=
% %             mnvl=mean(vls(indcheck))
% %           end
% %          
% %          indpass=find(abs(vls(indcheck)-mnvl)<sdpass&abs(vls(indcheck))>1.5);
% %         
% %          asympind=ind(indcheck(indpass));
% %          %THIS NEEDS TO BE MODIFIED TO RETURN ACSF RUNS PROPERLY
% %          if(~isempty(asympind))
% %             min_tm=floor(avls.tmvc(asympind(1),1));
% %             max_tm=floor(avls.tmvc(asympind(end),2));
% %             asympindall=find(avls.tmvc(:,1)>min_tm&avls.tmvc(:,1)<(max_tm+1))
% %          else
% %              asympindall=[];
% %          end
% %          
% %          asympvl=mean(aczvls(ind(indcheck(indpass))))
% %          muasympdist=asympvl-muzvls(ind(indcheck(indpass)));
% %          acasympdist=asympvl-aczvls(ind(indcheck(indpass)));
% %          if(asympvl<0)
% %              muasympdist=-muasympdist
% %              acasympdist=-acasympdist
% %          end
%   function [asympind,asympvl]=calcrevasympind(aczvls, muzvls,ind,bs)
%         
%          %start off with sdpass fixed to one
%          %i.e. calling an asymptote any set of consecutive days,
%          %at which <1sd from baseline.
%          sdpass=1;
%        
%           if ~isempty(bs.revasympvl)
%              mnvl=bs.revasympvl;
%              indcheck=1:length(aczvls(ind))
%              
%           else
% %             indcheck=find(abs(mxshift-vls)<sdnet)
% %             indcheckall=
%             mnvl=0;
%             indcheck=1:length(aczvls(ind))
%           end
%           
%           vls=aczvls(ind);
%          
%          indpass=find(abs(vls(indcheck)-mnvl)<sdpass);
%         
%          asympind=ind(indcheck(indpass));
%          asympvl=mean(aczvls(ind(indcheck(indpass))))
%          muasympdist=asympvl-muzvls(ind(indcheck(indpass)));
%          acasympdist=asympvl-aczvls(ind(indcheck(indpass)));
%          if(asympvl<0)
%              muasympdist=-muasympdist
%              acasympdist=-acasympdist
%          end
%          
%                  
%          