  %modified 4.11.10
  %purpose of this code is to go through sumbs,or phsumbs,and calculate
%summary output values, combvls, 
%in this version I take out all paired calculations.
function [shiftplot,combvls]=plothistoryscatter5(sumbs,ps)
DIFF=0
if(isfield(ps,'DIFF'))
    if(ps.DIFF)
        DIFF=1;
    end
end
if(isfield(ps,'USEPRE'))
    if(ps.USEPRE)
        USEPRE=1;
    else
        USEPRE=0;
    end
else
    USEPRE=0;
end

if(isfield(ps,'REV'))
    if(ps.REV)
        REV=1;
    else
        REV=0;
    end
else
    REV=0;
end

if(isfield(ps,'STIM'))
    if(ps.STIM)
        STIM=1;
    else
        STIM=0;
    end
else
    STIM=0;
end

NTANAL=ps.NTANAL

ctvl=1
 %1 is baseline %2 shift %3 rev % 4 restricted rev
    for colind=1:5
        for drxn=1:2
            for ntvl=NTANAL
            
                combvls{colind}{drxn}{ntvl}.offz=[]
                combvls{colind}{drxn}{ntvl}.mun=[];
                combvls{colind}{drxn}{ntvl}.acn=[];
                combvls{colind}{drxn}{ntvl}.acshift=[];
                combvls{colind}{drxn}{ntvl}.mushift=[];
                combvls{colind}{drxn}{ntvl}.acpreshift=[];
                combvls{colind}{drxn}{ntvl}.acpostshift=[];
                combvls{colind}{drxn}{ntvl}.pct=[];
                combvls{colind}{drxn}{ntvl}.bsnm=[];
                combvls{colind}{drxn}{ntvl}.shnum=[];
                combvls{colind}{drxn}{ntvl}.asymptms=[];
                combvls{colind}{drxn}{ntvl}.preinds=[];
                combvls{colind}{drxn}{ntvl}.pretms=[];
                combvls{colind}{drxn}{ntvl}.preindsh=[];
                combvls{colind}{drxn}{ntvl}.preindbs=[];
                combvls{colind}{drxn}{ntvl}.initshiftind=[];
                    combvls{colind}{drxn}{ntvl}.cvred=[];
                     combvls{colind}{drxn}{ntvl}.adjtimes=[];
                      combvls{colind}{drxn}{ntvl}.inds=[];
                        combvls{colind}{drxn}{ntvl}.acn=[];
                        combvls{colind}{drxn}{ntvl}.mun=[];
                
            end
        end
    end
for bsnm=1:length(sumbs)
    crbs=sumbs(bsnm);
    if(length(crbs.allnote)>1)
        ntind=crbs.allnote(1);
        ctrlind=crbs.allnote(2)
    else
        ntind=crbs.allnote(1);
        ctrlind=[];
    end
    shiftruns=crbs.shiftruns;
    
    revruns=crbs.revruns;
    numshiftruns=length(shiftruns);
    numrevruns=length(revruns);
    for ii=1:numshiftruns
        [EXCLUDE]=findexclude(bsnm,ii,ps);
        if(~EXCLUDE)
            asympruns=crbs.asympruns{ii};
        cnscasympruns=crbs.cnscasympruns{ii};
        if(~isempty(crbs.revasympruns))
            if(length(crbs.revasympruns)>=ii)
                revasympruns=crbs.revasympruns{ii};
                cnscrevasympruns=crbs.cnscrevasympruns{ii};
                revasympvl=crbs.revasympvl{ii};
            else
                revasympruns=[];
                cnscrevasympruns=[];
            end
        else
            revasympruns=[];
            cnscrevasympruns=[];
        end
        tms=crbs.adjshifttms{ii};
        shiftplot(bsnm).drxn{ii}=crbs.drxn{ii};
        %runtype is whether shift(1), bas(2), rev(3)
        for runtype=1:4
            if(runtype==1)
                inruns=crbs.shiftruns;
                
            elseif(runtype==2)
                inruns=crbs.sepbasruns;
            elseif(runtype==3)
                inruns=crbs.revruns;
            else
                
                inruns=crbs.asympruns;
                
            end
            if(ii<=length(inruns))
                     cr_runs{runtype}=inruns{ii};
                if(STIM)
                    if(runtype~=4)
                    tminds{runtype}=find(ismember(cr_runs{runtype},crbs.STANRUNS));
                     cr_runs{runtype}=cr_runs{runtype}(tminds{runtype});
                    else
                        tminds{runtype}=find(ismember(cr_runs{runtype},crbs.STANRUNS));
                        cr_runs{runtype}=cr_runs{runtype}(tminds{runtype});
      
                       
                        crtms=crbs.adjshifttms{ii}(tminds{1});
                        [preinds,pretms]=getpreruns(crtms,cr_runs);
                    end
                else
                    if(runtype~=4)
                    tminds{runtype}=find(ismember(cr_runs{runtype},crbs.mulist));
                     cr_runs{runtype}=cr_runs{runtype}(tminds{runtype});
                    else
                         tminds{runtype}=find(ismember(cr_runs{runtype},crbs.mulist));
                        cr_runs{runtype}=cr_runs{runtype}(tminds{runtype});
                        crtms=crbs.adjshifttms{ii}(tminds{1});
                        [preinds,pretms]=getpreruns(crtms,cr_runs);
                    end
                    
                end
                
                    
              tst=1  
            else
               cr_runs{runtype}=[]; 
            end
        end
    
        combins=[cr_runs{1}' cr_runs{2} ]
        shiftins=cr_runs{1};
        revinds=cr_runs{3}
        asympinds=cr_runs{4};
        if(STIM)
            aczcombvls{1}=crbs.acz
            acvls{1}=(crbs.acmean-crbs.initmean)./(crbs.sfact(1)/1000)
            muvls{1}=(crbs.mumean-crbs.initmean)./(crbs.sfact(1)/1000)
            acraw{1}=(crbs.acmean-crbs.initmean)./(crbs.sfact(1)/1000);
            muraw{1}=(crbs.mumean-crbs.initmean)./(crbs.sfact(1)/1000);
            muzvls{1}=crbs.muz
            cvred{1}=crbs.mu_cv./crbs.ac_cv;
            accombvls{1}=crbs.acmean;
           
        else
          
            for ntvl=ps.NTANAL
                if ntvl==1
                    nt=ntind;
                else
                    nt=ctrlind;
                end
                if(~isempty(crbs.acz(nt,:)))
                    aczprevls{ntvl}=crbs.acprez(nt,:);
                    crnt=find(crbs.allnote==nt);
                    acvls{ntvl}=(crbs.allmean(:,nt)-crbs.initmean(nt))./(crbs.sfact(crnt)/1000);
                    aczcombvls{ntvl}=crbs.acz(nt,:);
                    allraw{ntvl}=(crbs.allmean(:,nt)-crbs.initmean(nt))./(crbs.sfact(crnt)/1000);
                    
                    muzvls{ntvl}=crbs.muz(nt,:);
                    muvls{ntvl}=(crbs.allmean(:,nt)-crbs.initmean(nt))./(crbs.sfact(crnt)/1000);
                    aczpstvls{ntvl}=crbs.acpstz(nt,:);
                    cvred{ntvl}=calcpharmcv(crbs,nt);
                else
                    aczprevls{ntvl}=[];
                    aczcombvls{ntvl}=[];
                    acvls{ntvl}=[];
                    muzvls{ntvl}=[];
                    allraw{ntvl}=[];
                    muvls{ntvl}=[];
                    aczpstvls{ntvl}=[];
                    cvred{ntvl}=[];
                end
            
            end
        end
        %DIVIDE DATA INTO THREE GROUPS
        %BASIND(1)  SHIFTIND(2)  REVIND(3)  ASYMPIND(4)
        %pick the inds for baseline measurements, and shift measurements
        %and reverse to baseline, and reverse shift.
        analshiftinds=find(abs(aczcombvls{1}(shiftins))>1);
        allinds=find(abs(aczcombvls{1}(combins))>1);
        crtms=crbs.adjshifttms{ii};
        if(REV)
            shiftins=revinds
            analshiftinds=shiftinds;
            asympruns=revasympruns;
            cnscasympruns=cnscrevasympruns;
            crtms=crbs.revshifttms{ii}
        end
        
        analbasinds=setdiff(1:length(combins),allinds);
       
        crtms=crbs.adjshifttms{ii}
        %tminds{1} because this is shift run.
        mxtmind=find(crtms(tminds{1})<=ps.MX_TM);
        analshiftinds=intersect(mxtmind,analshiftinds);
        analshiftinds=shiftins(analshiftinds);
%         outinds{1}=combins(analbasinds);
        outinds{1}=cr_runs{2};
        %now remove asympruns from this analysis
        %and remove day max
        if(STIM)
            if(~isempty(asympruns))
                tmpind=find(ismember(asympruns,crbs.STANRUNS));
                asympruns=asympruns(tmpind);
            end
        end
        if(length(asympruns)>1)
          if(STIM)
              cnscruns=asympruns(2:end)
          else
              cnscruns=cnscasympruns;
          end
           if(~isempty(cnscruns))
                keepind=find(~ismember(analshiftinds,cnscruns))
           else
               keepind=1:length(analshiftinds);
           end
        end
           %throw out any runs after first asymp run.
           
           if(length(asympruns)>1)
                
               ntasympind=find(crbs.flrtmvec(analshiftinds)<=crbs.flrtmvec(asympruns(1)));
           end
               ntasympind=1:length(analshiftinds);
           
                combind=intersect(keepind,ntasympind);
           if(isempty(combind))
              combind=1; 
           end
           
%             keepind=analshiftinds;
           outinds{2}=analshiftinds; 
           allshiftinds=analshiftinds;
           
           tmindsout{2}=tminds{1};
           tmindsout{1}=tminds{2};
           tminds{1}=tmindsout{1};
           tminds{2}=tmindsout{2};
           allshift_tms=tmindsout{2};
           initshiftind=zeros(1,length(tminds{2}));
           initshiftind(combind)=1;
           %         
%         else
%             outinds{2}=analshiftinds;
%         end
        outinds{4}=asympinds;
        if(~isnan(preinds))
            outinds{5}=cr_runs{1}(preinds);
        else
            outinds{5}=[];
        end
        
        if(~isempty(cr_runs{3}))
%             analrevshiftinds=find(abs(aczvls(revind))>1);
%             analrevbasinds=setdiff(1:length(revind), analrevshiftinds);
%             analrevshiftinds=revind(analrevshiftinds);
%               if(length(revasympruns)>1)
                    if(~isempty(cnscrevasympruns))
                        keepind=find(~ismember(revinds,cnscrevasympruns));
                    else
                        keepind=1:length(revinds)
                    end
                    if(~isempty(revasympruns))
                        tmind=find(makerow((crbs.flrtmvec(revinds)<=crbs.flrtmvec(revasympruns(1)))));
                    else
                        tmind=1:length(revinds)
                    end
                        combind=intersect(keepind,tmind);
                    if(isempty(combind))
                        combind=1;
                    end
                    outinds{3}=revinds(combind);
                    tmpinds=find(abs(aczcombvls{1}(revinds(combind)))>1);
                    
                    tst=1;
%               else
%                   outinds{3}=revinds
%               end
        else
            outinds{3}=[];
            
        end

       if(~STIM)
        for runvl=1:5
            if(runvl~=3&runvl~=5)
                crtms=crbs.adjshifttms{ii}; 
            
            end
            for ntvl=ps.NTANAL
                if(~isempty(aczprevls{ntvl}))
                    if(ps.calcz)
                        %INSERT PRETMS AND PREINDS.
                        mulistinds=find(ismember(crbs.mulist,outinds{runvl}));
                        aclistinds=crbs.aclist(mulistinds,1);
                        aclistindspst=crbs.aclist(mulistinds,2);
                        shiftplot(bsnm).acdiff{ii}{runvl}{ntvl}=allraw{ntvl}(aclistinds);
                        if(runvl==2)
                            shiftplot(bsnm).initshiftind{ii}{runvl}{ntvl}=initshiftind
                        end
%                         shiftplot(bsnm).mun=crbs.mun(
                        shiftplot(bsnm).mun{ii}{runvl}{ntvl}=crbs.mun(ntvl,crbs.mulist(mulistinds));
                        shiftplot(bsnm).acn{ii}{runvl}{ntvl}=crbs.acn(ntvl,crbs.mulist(mulistinds));
                        shiftplot(bsnm).inds{ii}{runvl}{ntvl}=crbs.mulist(mulistinds);
                        shiftplot(bsnm).mudiff{ii}{runvl}{ntvl}=allraw{ntvl}(outinds{runvl});
                        shiftplot(bsnm).acdiffpst{ii}{runvl}{ntvl}=allraw{ntvl}(aclistindspst)
                        shiftplot(bsnm).acpreshift{ii}{runvl}{ntvl}=aczprevls{ntvl}(outinds{runvl});
                        shiftplot(bsnm).mushift{ii}{runvl}{ntvl}=muzvls{ntvl}(outinds{runvl});
                        shiftplot(bsnm).acpstshift{ii}{runvl}{ntvl}=aczpstvls{ntvl}(outinds{runvl});
                        
                        shiftplot(bsnm).accombshift{ii}{runvl}{ntvl}=aczcombvls{ntvl}(outinds{runvl});          
                        if(runvl<3)
                        shiftplot(bsnm).adjtms{ii}{runvl}{ntvl}=crtms(tminds{runvl});
                        end
                         shiftplot(bsnm).cvred{ii}{runvl}{ntvl}=cvred{ntvl}(outinds{runvl});
                        if(runvl==4)
                            
                            shiftplot(bsnm).asymptms{ii}{runvl}{ntvl}=(crbs.asympshifttms{ii}(tminds{4})-min(crbs.asympshifttms{ii})+1);
                            
                        end
                        if(runvl==5)
                            if(~isempty(preinds))
                            shiftplot(bsnm).preinds{ii}{runvl}{ntvl}=preinds;
                            shiftplot(bsnm).pretms{ii}{runvl}{ntvl}=pretms;
                            else
                                shiftplot(bsnm).preinds{ii}{runvl}{ntvl}=[];
                            shiftplot(bsnm).pretms{ii}{runvl}{ntvl}=[];
                            end    
                        end
                    %FOR LID RUNS.
                    %SET THE LID FLAG HERE.
                    else
                        mulistinds=find(ismember(crbs.mulist,outinds{runvl}));
                        aclistinds=crbs.aclist(mulistinds,1);
                        acpstinds=crbs.aclist(mulistinds,2);
%                         aclistinds=crbs.aclist(aclistinds,1);
                        shiftplot(bsnm).acpreshift{ii}{runvl}{ntvl}=acvls{ntvl}(aclistinds);
                        
                        shiftplot(bsnm).mushift{ii}{runvl}{ntvl}=muvls{ntvl}(outinds{runvl});
                        shiftplot(bsnm).acpstshift{ii}{runvl}{ntvl}=acvls{ntvl}(acpstinds);
                        shiftplot(bsnm).accombshift{ii}{runvl}{ntvl}=aczcombvls{ntvl}(outinds{runvl});          
                        if(runvl==2)
                            shiftplot(bsnm).initshiftind{ii}{runvl}{ntvl}=initshiftind
                        end
                        
                        
                        if(runvl<3)
                        shiftplot(bsnm).adjtms{ii}{runvl}{ntvl}=crtms(tminds{runvl});
                        end
                        shiftplot(bsnm).cvred{ii}{runvl}{ntvl}=cvred{ntvl}(outinds{runvl});
                        
                         if(runvl==4)
                            
                            shiftplot(bsnm).asymptms{ii}{runvl}{ntvl}=crbs.asympshifttms{ii}(tminds{4})-min(crbs.asympshifttms{ii})+1;
                           
                        end
                        if(runvl==5)
                            shiftplot(bsnm).preinds{ii}{runvl}{ntvl}=preinds;
                             shiftplot(bsnm).pretms{ii}{runvl}{ntvl}=pretms;
                        end
                            
                    end
                    
                    
%                 shiftplot(bsnm).mushift{ii}{runvl}{ntvl}=muzvls{ntvl}(outinds{runvl});
                else
                    shiftplot(bsnm).acpreshift{ii}{runvl}{ntvl}=NaN;
                    shiftplot(bsnm).mushift{ii}{runvl}{ntvl}=NaN;
                     shiftplot(bsnm).cvred{ii}{runvl}{ntvl}=NaN;
                    shiftplot(bsnm).acpstshift{ii}{runvl}{ntvl}=NaN;
                    shiftplot(bsnm).accombshift{ii}{runvl}{ntvl}=NaN;          
                    shiftplot(bsnm).adjtms{ii}{runvl}{ntvl}=NaN;
%                     shiftplot(bsnm).pretms{ii}{runvl}{ntvl}=NaN;
%                     
                end
            tst=1;
            end
        end
       else
          for runvl=1:5
               if(runvl~=3&runvl~=5)
                crtms=crbs.adjshifttms{ii}; 
            
            end
            for ntvl=ps.NTANAL
                if(~isempty(aczcombvls{ntvl}))
                    if(ps.calcz)
                        shiftplot(bsnm).mushift{ii}{runvl}{ntvl}=muzvls{ntvl}(outinds{runvl});
                         shiftplot(bsnm).acdiff{ii}{runvl}{ntvl}=acraw{ntvl}(outinds{runvl});
                        shiftplot(bsnm).mudiff{ii}{runvl}{ntvl}=muraw{ntvl}(outinds{runvl});
                         if(runvl==2)
                            shiftplot(bsnm).initshiftind{ii}{runvl}{ntvl}=initshiftind
                        end
                        
                        
                        if(runvl~=3&runvl~=5)
                           
                        shiftplot(bsnm).adjtms{ii}{runvl}{ntvl}=crbs.adjshifttms{ii}(tminds{runvl});
                        end
                        shiftplot(bsnm).accombshift{ii}{runvl}{ntvl}=aczcombvls{ntvl}(outinds{runvl});  
                        shiftplot(bsnm).cvred{ii}{runvl}{ntvl}=cvred{ntvl}(outinds{runvl});
                         if(runvl==4)
                            
                             shiftplot(bsnm).asymptms{ii}{runvl}{ntvl}=crbs.asympshifttms{ii}(tminds{4})-min(crbs.asympshifttms{ii})+1;
                            
                         end


                         if(runvl==5)
                            shiftplot(bsnm).preinds{ii}{runvl}{ntvl}=preinds;
                            shiftplot(bsnm).pretms{ii}{runvl}{ntvl}=pretms;
                         end
                    else
                        shiftplot(bsnm).mushift{ii}{runvl}{ntvl}=muvls{ntvl}(outinds{runvl});
                         if(runvl~=3&runvl~=5)
                        shiftplot(bsnm).adjtms{ii}{runvl}{ntvl}=crbs.adjshifttms{ii}(tminds{runvl});
                         end
                        shiftplot(bsnm).accombshift{ii}{runvl}{ntvl}=acvls{ntvl}(outinds{runvl}); 
                        shiftplot(bsnm).acdiff{ii}{runvl}{ntvl}=acraw{ntvl}(outinds{runvl});
                        shiftplot(bsnm).mudiff{ii}{runvl}{ntvl}=muraw{ntvl}(outinds{runvl});
                        shiftplot(bsnm).cvred{ii}{runvl}{ntvl}=cvred{ntvl}(outinds{runvl});
                     if(runvl==2)
                            shiftplot(bsnm).initshiftind{ii}{runvl}{ntvl}=initshiftind
                        end
                        
                        if(runvl==4) 
                             shiftplot(bsnm).asymptms{ii}{runvl}{ntvl}=crbs.asympshifttms{ii}(tminds{4})-min(crbs.asympshifttms{ii})+1;
                        end
                    
                     if(runvl==5)
                            shiftplot(bsnm).preinds{ii}{runvl}{ntvl}=preinds;
                            shiftplot(bsnm).pretms{ii}{runvl}{ntvl}=pretms;
                         end
                    end
                else
                   
                    shiftplot(bsnm).mushift{ii}{runvl}{ntvl}=NaN;
                    shiftplot(bsnm).accombshift{ii}{runvl}{ntvl}=NaN;         
                    shiftplot(bsnm).cvred{ii}{runvl}{ntvl}=NaN;
                    shiftplot(bsnm).adjtms{ii}{runvl}{ntvl}=crtms{ii}(tminds);
                end
            end
          end 
       end
       end
    end 
    ps.bsnm=bsnm;
%     ps.DIFF=0;
    [combvls]=meanhorizarrow(shiftplot(bsnm),combvls,ps);
    ctvl=ctvl+1;
    
end
% linkaxes(ax);

    function [preinds,pretms]=getpreruns(tminds,cr_runs)
    preinds=find(~ismember(cr_runs{1},cr_runs{4}));
    if(~isempty(cr_runs{4}))
        asympinitind=find(ismember(cr_runs{1},cr_runs{4}));
        asymptm=tminds(asympinitind(1));
        if(~isempty(preinds))
            pretms=tminds(preinds)-asymptm;
        else
            pretms=[]
        end
        outind=find(pretms<=0);
        pretms=pretms(outind);
        preinds=preinds(outind);
    else
        pretms=NaN;
        preinds=NaN;
    end
    

    function [cvredout]=calcpharmcv(crbs,nt)
        mulistind=find(crbs.mulist>0);
        for ii=1:length(mulistind)
            crind=mulistind(ii);
           muind=crbs.mulist(crind);
           acind=crbs.aclist(crind,1);
           cvredout(muind)=crbs.ac_cv(nt,muind)./crbs.ac_cv(nt,acind);
            
        end
    function [EXCLUDE]=findexclude(bsnm,insh,ps)
        excludeind=0;
        EXCLUDE=0;
        ct=0;
        if(isfield(ps,'excludeind'))
            if(ps.excludeind)
                excludeind=1;
            end
        end
        if(excludeind)
            for ii=1:length(ps.excludebs)
                crbs=ps.excludebs{ii};
                crsh=ps.excludeshnum{ii};
                
                if(bsnm==crbs&insh==insh)
                    ct=ct+1;
                else
                    ct=0;
                end
            end
        end
        if(ct)
            EXCLUDE=1;
        end
        
function [combvls]=meanhorizarrow(shiftplot,combvls,ps)
% axes(ax);
bsnm=ps.bsnm;
hold on;
numdrxn=length(shiftplot.drxn);  
% plot([-4 4],[ctvl-.3 ctvl-.3],'Color',[0.6 0.6 0.6])
for ii=1:numdrxn
    if(ps.USEPRE)
       if(~ps.DIFF)
            crac=shiftplot.acpreshift{ii}
       else
           crac=shiftplot.acdiff{ii}
       end
    else
       if(~ps.DIFF)
            crac=shiftplot.accombshift{ii}
       else
           crac=shiftplot.acdiff{ii}
       end
            
    end
      if(~ps.DIFF)
        crmu=shiftplot.mushift{ii}
      else
          crmu=shiftplot.mudiff{ii}
      end
       asymptms=shiftplot.asymptms{ii};
       if(ps.STIM==0)
            cracpst=shiftplot.acpstshift{ii}
            cr_cvred=shiftplot.cvred{ii};
       else
           cr_cvred=shiftplot.cvred{ii};
       end
       for typevl=1:5
        for ntvl=ps.NTANAL
            drxn=shiftplot.drxn{ii}
%             acn=shiftplot.acn{ii}{typevl}{ntvl};
%             mun=shiftplot.mun{ii}{typevl}{ntvl};
%             inds=shiftplot.inds{ii}{typevl}{ntvl};
           
            if(drxn=='up')
                drxnvl=1;
            else
                drxnvl=2;
            end
            numvls=length(crac{typevl});
           if(ntvl<=length((shiftplot.preinds{ii}{5})))
                preinds=shiftplot.preinds{ii}{5}{ntvl}
                pretms=shiftplot.pretms{ii}{5}{ntvl}
           else
               preinds=[];
               pretms=[];
           end
             adjtms=shiftplot.adjtms{ii}{2}{ntvl}
            
            if(numvls)
                if(~isnan(crac{typevl}{ntvl}))
                    y1=crac{typevl}{ntvl};
                    
                    if(typevl==2)
                        initshiftind=shiftplot.initshiftind{ii}{2}{ntvl}
                    end
                    y2=crmu{typevl}{ntvl};
                    if(typevl==4)
                        crasymptms=asymptms{typevl}{ntvl}
                    end
                    
                        cr_cvredout=cr_cvred{typevl}{ntvl}
                    
                    diff=y2-y1;
               
                    if(drxn=='up')
                        combdiff=diff
                        
                    else
                        combdiff=diff;
                        
                    end
                    pct=diff./y1
                    nonaninds=find(~isnan(combdiff))
                else
                    nonaninds=0;
                end
                if(typevl==5)
                    nonaninds=1:length(y1)
                end
                
                    %TARGET NOTE
                    if(nonaninds)
                        combvls{typevl}{drxnvl}{ntvl}.offz=[makerow(combvls{typevl}{drxnvl}{ntvl}.offz) makerow(combdiff(nonaninds))]
                        combvls{typevl}{drxnvl}{ntvl}.acshift=[makerow(combvls{typevl}{drxnvl}{ntvl}.acshift) makerow(y1(nonaninds))]
                        combvls{typevl}{drxnvl}{ntvl}.mushift=[makerow(combvls{typevl}{drxnvl}{ntvl}.mushift) makerow(y2(nonaninds))];
                        combvls{typevl}{drxnvl}{ntvl}.pct=[makerow(combvls{typevl}{drxnvl}{ntvl}.pct) makerow(pct(nonaninds))];
                        combvls{typevl}{drxnvl}{ntvl}.bsnm=[combvls{typevl}{drxnvl}{ntvl}.bsnm bsnm*ones(1,length(nonaninds))];
                        combvls{typevl}{drxnvl}{ntvl}.shnum=[combvls{typevl}{drxnvl}{ntvl}.shnum ii*ones(1,length(nonaninds))];
                        combvls{typevl}{drxnvl}{ntvl}.cvred=[combvls{typevl}{drxnvl}{ntvl}.cvred cr_cvredout(nonaninds)];
                        if(typevl<5)
%                         combvls{typevl}{drxnvl}{ntvl}.inds=[combvls{typevl}{drxnvl}{ntvl}.inds inds(nonaninds)];
%                         combvls{typevl}{drxnvl}{ntvl}.acn=[combvls{typevl}{drxnvl}{ntvl}.acn acn(nonaninds)];
%                         combvls{typevl}{drxnvl}{ntvl}.mun=[combvls{typevl}{drxnvl}{ntvl}.mun mun(nonaninds)];
                        end
                        if(typevl==2)
                            combvls{typevl}{drxnvl}{ntvl}.initshiftind=[combvls{typevl}{drxnvl}{ntvl}.initshiftind initshiftind(nonaninds)];
                        end
                        if(typevl==4)
                            combvls{typevl}{drxnvl}{ntvl}.asymptms=[combvls{typevl}{drxnvl}{ntvl}.asymptms makerow(crasymptms(nonaninds))]
                            combvls{typevl}{drxnvl}{ntvl}.preinds=[combvls{typevl}{drxnvl}{ntvl}.preinds makerow(preinds)];
                            combvls{typevl}{drxnvl}{ntvl}.pretms=[combvls{typevl}{drxnvl}{ntvl}.pretms makerow(pretms)]
                            combvls{typevl}{drxnvl}{ntvl}.preindsh=[combvls{typevl}{drxnvl}{ntvl}.preindsh ii*ones(1,length(pretms))]
                            combvls{typevl}{drxnvl}{ntvl}.preindbs=[combvls{typevl}{drxnvl}{ntvl}.preindbs bsnm*ones(1,length(pretms))]
                        end
                            if(ps.STIM)
                            combvls{typevl}{drxnvl}{ntvl}.adjtimes=[makerow(combvls{typevl}{drxnvl}{ntvl}.adjtimes) makerow(adjtms(nonaninds))] 
                            combvls{typevl}{drxnvl}{ntvl}.acpreshift=[makerow(combvls{typevl}{drxnvl}{ntvl}.acpreshift) makerow(y1(nonaninds))]
                            else
                            if(typevl<5)
                                combvls{typevl}{drxnvl}{ntvl}.adjtimes=[makerow(combvls{typevl}{drxnvl}{ntvl}.adjtimes) makerow(adjtms(nonaninds))] 
                                combvls{typevl}{drxnvl}{ntvl}.acpreshift=[makerow(combvls{typevl}{drxnvl}{ntvl}.acpreshift) makerow(y1(nonaninds))]
                                combvls{typevl}{drxnvl}{ntvl}.acpostshift=[makerow(combvls{typevl}{drxnvl}{ntvl}.acpostshift) makerow(cracpst{typevl}{ntvl}(nonaninds))];
                            end
                            end
                    else
                        combvls{typevl}{drxnvl}{ntvl}.offz=[combvls{typevl}{drxnvl}{ntvl}.offz NaN]
                       
                         if(typevl==2)
                            combvls{typevl}{drxnvl}{ntvl}.initshiftind=[combvls{typevl}{drxnvl}{ntvl}.initshiftind NaN];
                        end
                        
                        if(typevl<5)
                        combvls{typevl}{drxnvl}{ntvl}.acshift=[combvls{typevl}{drxnvl}{ntvl}.acshift NaN]
%                         
                        combvls{typevl}{drxnvl}{ntvl}.mushift=[combvls{typevl}{drxnvl}{ntvl}.mushift NaN];
                       end
                         combvls{typevl}{drxnvl}{ntvl}.pct=[combvls{typevl}{drxnvl}{ntvl}.pct NaN];
                        combvls{typevl}{drxnvl}{ntvl}.bsnm=[combvls{typevl}{drxnvl}{ntvl}.bsnm bsnm];
                        combvls{typevl}{drxnvl}{ntvl}.shnum=[combvls{typevl}{drxnvl}{ntvl}.shnum ii];
                        combvls{typevl}{drxnvl}{ntvl}.cvred=[combvls{typevl}{drxnvl}{ntvl}.cvred NaN];
                        if(typevl==4)
                            combvls{typevl}{drxnvl}{ntvl}.asymptms=[combvls{typevl}{drxnvl}{ntvl}.asymptms NaN]
%                             combvls{typevl}{drxnvl}{ntvl}.preinds=[combvls{typevl}{drxnvl}{ntvl}.preinds NaN];
%                             combvls{typevl}{drxnvl}{ntvl}.pretms=[combvls{typevl}{drxnvl}{ntvl}.pretms NaN]
%                             combvls{typevl}{drxnvl}{ntvl}.preindsh=[combvls{typevl}{drxnvl}{ntvl}.preindsh NaN]
%                             combvls{typevl}{drxnvl}{ntvl}.preindbs=[combvls{typevl}{drxnvl}{ntvl}.preindbs NaN]
                        end     
                        combvls{typevl}{drxnvl}{ntvl}.adjtimes=[makerow(combvls{typevl}{drxnvl}{ntvl}.adjtimes) NaN] 
                        if(ps.STIM)
                            
                            combvls{typevl}{drxnvl}{ntvl}.acpreshift=[makerow(combvls{typevl}{drxnvl}{ntvl}.acpreshift) NaN]
                        else
                            combvls{typevl}{drxnvl}{ntvl}.acpreshift=[combvls{typevl}{drxnvl}{ntvl}.acpreshift NaN]
                            combvls{typevl}{drxnvl}{ntvl}.acpostshift=[combvls{typevl}{drxnvl}{ntvl}.acpostshift NaN];
                        end
                     end  
                        %CONTROL
%                     else
%                         combvls{1}{drxnvl}.acshiftctrl=[combvls{1}{drxnvl}.acshiftctrl y1ctrl(nonaninds)]
%                         combvls{1}{drxnvl}.mushiftctrl=[combvls{1}{drxnvl}.mushiftctrl y2ctrl(nonaninds)]
%                     else
%                        combvls{1}{drxnvl}.acshiftctrl=[combvls{1}{drxnvl}.acshiftctrl nan(1,length(nonaninds))]
%                         combvls{1}{drxnvl}.mushiftctrl=[combvls{1}{drxnvl}.mushiftctrl nan(1,length(nonaninds))] 
%                     end
                       
               
        end
        end
    end
end


