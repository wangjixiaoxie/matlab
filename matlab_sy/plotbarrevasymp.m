%written 4.1.10 to plot group data histograms
%input combvls takes the following form.
%combvls{1} - baseline runs
%combvls{2} - shift runs

% %REWRITING FROM SCRATCH - rather than trying to modify existing
% shiftplot, combvls code
function [rs_stim,rs_pharm,plotvls_stim,plotvls_pharm]=plotbarrevasymp(sumbs,phsumbs)   
    %first do analysis for stim.
    %need structure which lists tmvalues,acvalues,muvalues
    %bs, sh, drxn
%     %asympvalues
    ps.STIM=1;
     ps.tmwins{1}=[0];
    ps.tmwins{2}=[2:3]
    [rs_stim]=getvals(sumbs,ps);
    
     [plotvls_stim]=calcplotvls(sumbs,rs_stim,ps);
    
    ps.STIM=0;
    [rs_pharm]=getvals(phsumbs,ps);
    ps.tmwins{1}=[0];
    ps.tmwins{2}=[2:3]
%    ps.tmwins{3}=[3:4]; 
    
    [plotvls_pharm]=calcplotvls(phsumbs,rs_pharm,ps);
    
 function [plotvls]=calcplotvls(sumbs,rs,ps)
       %get unique bs
       plotvls.mndiffhz=[];
%        plotvls.meanmu=[];
       plotvls.crbs=[];
       plotvls.crsh=[];
       
       unqbs=unique(rs.bs);
       for ii=1:length(unqbs)
          crbs=unqbs(ii);
          crbsinds=find(rs.bs==crbs);
          unqsh=unique(rs.sh(crbsinds));
          for jj=1:length(unqsh)
              crsh=unqsh(jj);
              crinds=find(rs.bs==crbs&rs.sh==crsh);
              revtms=rs.revtms(crinds);
              asympflag=rs.asympflag(crinds);
              asympind=find(asympflag==1);
              if(~isempty(asympind))
                asymptms=revtms(asympind)-revtms(asympind(1));
                ps.cr_asympinds=rs.revruns((crinds(asympind)));
                TM_MATCH=testtmmatch(asymptms,ps);
                ps.crbs=crbs;
                ps.crsh=crsh;
              else
                  TM_MATCH=0;
              end
              if(TM_MATCH)
                  [mndiffhz]=getplotvls(sumbs,rs,asympind,asymptms,ps);
                  plotvls.mndiffhz=[plotvls.mndiffhz;mndiffhz];
%                   plotvls.meanmu=[plotvls.meanmu;meanmu];
                  plotvls.crbs=[plotvls.crbs crbs];
                  plotvls.crsh=[plotvls.crsh crsh];
              end
          end
       end

function [mndiffhz]=getplotvls(sumbs,rs,asympind,asymptms,ps)
   for crtmind=1:length(ps.tmwins)
      	crinds=ps.cr_asympinds;
        crtms=ps.tmwins{crtmind}
        [vls,matchinds]=intersect(asymptms,crtms);
%         inds=revinds(asympind(matchinds));
        inds=crinds(matchinds);
        %THIS NEEDS TO BE WORKED OUT.
        
        [mndiffhz(crtmind)]=getacmu(sumbs,ps,inds);
        
   end
function [mndiffhz]=getacmu(sumbs,ps,inds)
    if(sumbs(ps.crbs).drxn{ps.crsh}=='do')
        FACT=1;
    else
        FACT=-1;
    end
    if(ps.STIM)
         crbs=sumbs(ps.crbs);
        ntind=crbs.ntind;
        sfact=crbs.sfact(1);
         rawdiffhz=FACT*(sumbs(ps.crbs).acmean(inds)-sumbs(ps.crbs).mumean(inds))./(sfact/1000);
         mndiffhz=nanmean(rawdiffhz);   
        
      tst=1;  
    else
        crbs=sumbs(ps.crbs);
        ntind=crbs.ntind;
        sfact=crbs.sfact(1);
        indloc=find(ismember(crbs.mulist,inds));
        acloc=crbs.aclist(indloc,1);
        acvls=crbs.allmean(acloc,ntind);
        muvl=crbs.allmean(inds,ntind);
        rawdiffhz=FACT*(acvls-muvl)./(sfact/1000);
       mndiffhz=nanmean(rawdiffhz);
    end
    

    function [TM_MATCH]=testtmmatch(asymptms,ps)
       
           for crtmind=1:length(ps.tmwins)
                   
               crtms=ps.tmwins{crtmind}
               [vls,matchinds]=intersect(asymptms,crtms);
               if(~isempty(matchinds))
                   match(crtmind)=1;
               else
                   match(crtmind)=NaN;
               end
           end
           
           notnaind=find(isnan(match));
           if(isempty(notnaind))
               TM_MATCH=1;
           else
               TM_MATCH=0;
           end
   
       
   function [rs]=getvals(sumbs,ps)
       rs.bs=[];
       rs.sh=[];
       rs.revtms=[];
       rs.asympflag=[];
       rs.revruns=[];
       for ii=1:length(sumbs)
           crbs=sumbs(ii);
          for sh=1:length(crbs.revruns)
              asympflag=[];
            revruns=crbs.revruns{sh};
            allasympruns=crbs.revasympruns{sh};
            asymprevind=find(ismember(revruns,allasympruns));
            asympflag(1:length(revruns))=0;
            asympflag(asymprevind)=1;
            if(~isempty(crbs.revshifttms))
                if(sh<=length(crbs.revshifttms))
                    revtms=crbs.revshifttms{sh}
                    rs.revruns=[rs.revruns makerow(revruns)];
                    rs.asympflag=[rs.asympflag makerow(asympflag)];
                    rs.revtms=[rs.revtms makerow(revtms)];
                    rs.bs=[rs.bs ii*ones(1,length(revtms))];
                    rs.sh=[rs.sh sh*ones(1,length(revtms))];
          
                end
            end
          end
       end
            
       
       
 
