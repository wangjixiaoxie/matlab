  %modified 4.11.10
  %purpose of this code is to go through sumbs,or phsumbs,and calculate
%summary output values, combvls, 
%in this version I take out all paired calculations.
function [shiftplot,combvls]=plothistoryscatter4(sumbs,ps)
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
    for colind=1:4
        for drxn=1:2
            for ntvl=NTANAL
            
                combvls{colind}{drxn}{ntvl}.offz=[]
                combvls{colind}{drxn}{ntvl}.acshift=[];
                combvls{colind}{drxn}{ntvl}.mushift=[];
                combvls{colind}{drxn}{ntvl}.acpreshift=[];
                combvls{colind}{drxn}{ntvl}.acpostshift=[];
                combvls{colind}{drxn}{ntvl}.pct=[];
                combvls{colind}{drxn}{ntvl}.bsnm=[];
                combvls{colind}{drxn}{ntvl}.shnum=[];
                if(STIM)
                    combvls{colind}{drxn}{ntvl}.cvred=[];
                end
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
        for runtype=1:3
            if(runtype==1)
                inruns=crbs.shiftruns;
                
            elseif(runtype==2)
                inruns=crbs.sepbasruns;
            else
                inruns=crbs.revruns;
            end
            if(ii<=length(inruns))
                if(runtype==2)
                    cr_runs{runtype}=inruns{ii};
                else
                     cr_runs{runtype}=inruns{ii};
                end
                if(STIM)
                    tmpind=find(ismember(cr_runs{runtype},crbs.STANRUNS));
                    cr_runs{runtype}=cr_runs{runtype}(tmpind);
                end
            else
               cr_runs{runtype}=[]; 
            end
        end
    
        combins=[cr_runs{1}' cr_runs{2} ]
        shiftins=cr_runs{1};
        revinds=cr_runs{3}
        if(STIM)
            aczcombvls{1}=crbs.acz
            acvls{1}=(crbs.acmean-crbs.initmean)./(crbs.sfact(1)/1000)
            muvls{1}=(crbs.mumean-crbs.initmean)./(crbs.sfact(1)/1000)
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
                    muzvls{ntvl}=crbs.muz(nt,:);
                    muvls{ntvl}=(crbs.allmean(:,nt)-crbs.initmean(nt))./(crbs.sfact(crnt)/1000);
                    aczpstvls{ntvl}=crbs.acpstz(nt,:);
                else
                    aczprevls{ntvl}=[];
                    aczcombvls{ntvl}=[];
                    acvls{ntvl}=[];
                    muzvls{ntvl}=[];
                    muvls{ntvl}=[];
                    aczpstvls{ntvl}=[];
                end
            
            end
        end
        %DIVIDE DATA INTO THREE GROUPS
        %BASIND(1)  SHIFTIND(2)  REVIND(3)
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
        tminds=find(crtms<=ps.MX_TM);
        analshiftinds=intersect(tminds,analshiftinds);
        analshiftinds=shiftins(analshiftinds);
%         outinds{1}=combins(analbasinds);
        outinds{1}=cr_runs{2};
        %now remove asympruns from this analysis
        %and remove day max
        
%         if(length(asympruns)>1)
           %throw out consecutive runs
           if(~isempty(cnscasympruns))
                keepind=find(~ismember(analshiftinds,cnscasympruns))
           else
               keepind=1:length(analshiftinds);
           end
           %throw out any runs after first asymp run.
           if(~isempty(asympruns))
                tmind=find(crbs.flrtmvec(analshiftinds)<=crbs.flrtmvec(asympruns(1)));
           else
               tmind=1:length(analshiftinds);
           end
                combind=intersect(keepind,tmind);
           if(isempty(combind))
              combind=1; 
           end
           
%             keepind=analshiftinds;
           outinds{2}=analshiftinds(combind); 
%         else
%             outinds{2}=analshiftinds;
%         end
        
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
                    outinds{4}=revinds(combind(tmpinds));
                    tst=1;
%               else
%                   outinds{3}=revinds
%               end
        else
            outinds{3}=[];
            outinds{4}=[]
        end

       if(~STIM)
        for runvl=1:4
            for ntvl=ps.NTANAL
                if(~isempty(aczprevls{ntvl}))
                    if(ps.calcz)
                        shiftplot(bsnm).acpreshift{ii}{runvl}{ntvl}=aczprevls{ntvl}(outinds{runvl});
                        shiftplot(bsnm).mushift{ii}{runvl}{ntvl}=muzvls{ntvl}(outinds{runvl});
                        shiftplot(bsnm).acpstshift{ii}{runvl}{ntvl}=aczpstvls{ntvl}(outinds{runvl});
                        shiftplot(bsnm).accombshift{ii}{runvl}{ntvl}=aczcombvls{ntvl}(outinds{runvl});          
                        shiftplot(bsnm).adjtms{ii}{runvl}{ntvl}=crbs.adjshifttms{ii}(tminds);
                    %FOR LID RUNS.
                    %SET THE LID FLAG HERE.
                    else
                        mulistinds=find(ismember(crbs.mulist,outinds{runvl}));
                        aclistinds=crbs.aclist(mulistinds,1);
%                         aclistinds=crbs.aclist(aclistinds,1);
                        shiftplot(bsnm).acpreshift{ii}{runvl}{ntvl}=acvls{ntvl}(aclistinds);
                        shiftplot(bsnm).mushift{ii}{runvl}{ntvl}=muvls{ntvl}(outinds{runvl});
                        shiftplot(bsnm).acpstshift{ii}{runvl}{ntvl}=aczpstvls{ntvl}(outinds{runvl});
                        shiftplot(bsnm).accombshift{ii}{runvl}{ntvl}=aczcombvls{ntvl}(outinds{runvl});          
                        shiftplot(bsnm).adjtms{ii}{runvl}{ntvl}=crbs.adjshifttms{ii}(tminds);
                        
                        
                    end
                    
                    
%                 shiftplot(bsnm).mushift{ii}{runvl}{ntvl}=muzvls{ntvl}(outinds{runvl});
                else
                    shiftplot(bsnm).acpreshift{ii}{runvl}{ntvl}=NaN;
                    shiftplot(bsnm).mushift{ii}{runvl}{ntvl}=NaN;
                    shiftplot(bsnm).acpstshift{ii}{runvl}{ntvl}=NaN;
                    shiftplot(bsnm).accombshift{ii}{runvl}{ntvl}=NaN;          
                    shiftplot(bsnm).adjtms{ii}{runvl}{ntvl}=crbs.adjshifttms{ii}(tminds);
                end
            end
        end
       else
          for runvl=1:4
            for ntvl=ps.NTANAL
                if(~isempty(aczcombvls{ntvl}))
                    if(ps.calcz)
                        shiftplot(bsnm).mushift{ii}{runvl}{ntvl}=muzvls{ntvl}(outinds{runvl});
                   
                        shiftplot(bsnm).accombshift{ii}{runvl}{ntvl}=aczcombvls{ntvl}(outinds{runvl});  
                        shiftplot(bsnm).cvred{ii}{runvl}{ntvl}=cvred{ntvl}(outinds{runvl});
                    else
                        shiftplot(bsnm).mushift{ii}{runvl}{ntvl}=muvls{ntvl}(outinds{runvl});
                   
                        shiftplot(bsnm).accombshift{ii}{runvl}{ntvl}=acvls{ntvl}(outinds{runvl});  
                        shiftplot(bsnm).cvred{ii}{runvl}{ntvl}=cvred{ntvl}(outinds{runvl});
                    end
                else
                   
                    shiftplot(bsnm).mushift{ii}{runvl}{ntvl}=NaN;
                    shiftplot(bsnm).accombshift{ii}{runvl}{ntvl}=NaN;         
                    shiftplot(bsnm).cvred{ii}{runvl}{ntvl}=NaN;
                end
            end
          end 
       end
    end 
    ps.bsnm=bsnm;
    [combvls]=meanhorizarrow(shiftplot(bsnm),combvls,ps);
    ctvl=ctvl+1;
    
end
% linkaxes(ax);


function [combvls]=meanhorizarrow(shiftplot,combvls,ps)
% axes(ax);
bsnm=ps.bsnm;
hold on;
numdrxn=length(shiftplot.drxn);  
% plot([-4 4],[ctvl-.3 ctvl-.3],'Color',[0.6 0.6 0.6])
for ii=1:numdrxn
    if(ps.USEPRE)
       crac=shiftplot.acpreshift{ii}
       
       
    else
        crac=shiftplot.accombshift{ii}
    end
       crmu=shiftplot.mushift{ii} 
       
       if(ps.STIM==0)
            cracpst=shiftplot.acpstshift{ii}
       else
           cr_cvred=shiftplot.cvred{ii};
       end
       for typevl=1:4
        for ntvl=ps.NTANAL
            drxn=shiftplot.drxn{ii};   
            if(drxn=='up')
                drxnvl=1;
            else
                drxnvl=2;
            end
            numvls=length(crac{typevl});
            if(numvls)
                if(~isnan(crac{typevl}{ntvl}))
                    y1=crac{typevl}{ntvl};
                    y2=crmu{typevl}{ntvl};
                    if(ps.STIM)
                        cr_cvredout=cr_cvred{typevl}{ntvl}
                    end
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
                
                
                    %TARGET NOTE
                    if(nonaninds)
                        combvls{typevl}{drxnvl}{ntvl}.offz=[makerow(combvls{typevl}{drxnvl}{ntvl}.offz) makerow(combdiff(nonaninds))]
                        combvls{typevl}{drxnvl}{ntvl}.acshift=[makerow(combvls{typevl}{drxnvl}{ntvl}.acshift) makerow(y1(nonaninds))]
                        combvls{typevl}{drxnvl}{ntvl}.mushift=[makerow(combvls{typevl}{drxnvl}{ntvl}.mushift) makerow(y2(nonaninds))];
                        combvls{typevl}{drxnvl}{ntvl}.pct=[makerow(combvls{typevl}{drxnvl}{ntvl}.pct) makerow(pct(nonaninds))];
                        combvls{typevl}{drxnvl}{ntvl}.bsnm=[combvls{typevl}{drxnvl}{ntvl}.bsnm bsnm*ones(1,length(nonaninds))];
                        combvls{typevl}{drxnvl}{ntvl}.shnum=[combvls{typevl}{drxnvl}{ntvl}.shnum ii*ones(1,length(nonaninds))];
                        if(ps.STIM)
                            combvls{typevl}{drxnvl}{ntvl}.cvred=[combvls{typevl}{drxnvl}{ntvl}.cvred cr_cvredout(nonaninds)];
                        else
                            combvls{typevl}{drxnvl}{ntvl}.acpreshift=[makerow(combvls{typevl}{drxnvl}{ntvl}.acpreshift) makerow(y1(nonaninds))]
                            combvls{typevl}{drxnvl}{ntvl}.acpostshift=[makerow(combvls{typevl}{drxnvl}{ntvl}.acpostshift) makerow(cracpst{typevl}{ntvl}(nonaninds))];
                        end
                    else
                        combvls{typevl}{drxnvl}{ntvl}.offz=[combvls{typevl}{drxnvl}{ntvl}.offz NaN]
                        combvls{typevl}{drxnvl}{ntvl}.acshift=[combvls{typevl}{drxnvl}{ntvl}.acshift NaN]
                        
                        combvls{typevl}{drxnvl}{ntvl}.mushift=[combvls{typevl}{drxnvl}{ntvl}.mushift NaN];
                         
                         combvls{typevl}{drxnvl}{ntvl}.pct=[combvls{typevl}{drxnvl}{ntvl}.pct NaN];
                        combvls{typevl}{drxnvl}{ntvl}.bsnm=[combvls{typevl}{drxnvl}{ntvl}.bsnm bsnm];
                        combvls{typevl}{drxnvl}{ntvl}.shnum=[combvls{typevl}{drxnvl}{ntvl}.shnum ii];
                        if(ps.STIM)
                            combvls{typevl}{drxnvl}{ntvl}.cvred=[combvls{typevl}{drxnvl}{ntvl}.cvred NaN];
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


