%MAJOR Modifications and updating 5.19.09


%OUTPUTS
%sumplotout, sumplotrevout

%INPUTS
%range is two component cell array, e.g. {1:5 1:5}
%minshift is two component vector [1.5 1.5]
%modstruct useful for exclusion
%e.g. bk63w43
%modstruct.bsnum=9;
%modstruct.shftnum=1;
% % modstruct.runs=20;
% 
% modstruct.bsnum=9;
% modstruct.shftnum=1;
% modstruct.runs=20;
% 
% %for cell array indexing
%1 is target
%2 is control
%3 is baseline


function [sumplotout,sumplotrevout]=selectinbasruns4(sumbs,init_typ, range, minshift,modstruct,min_notes)

%init_typ 1 means days since wn was turned on.
%init_typ 2 means days since crossing threshold

%range is values, for example 1:5, includes first five days since wn was
%turned on.

%sumtype 1 is for sumplotout, sumplotrevout is sumtype 2
for sumtype=1
% for sumtype=1:2
    clear sumplot sumplotrevout
    
    datout=initfields();
    for bsnum=1:length(sumbs)
         crbs=sumbs(bsnum)
         if(exist('modstruct'))
            [initind, initday, drxnvec,basind, shftnum,asympvl]=getshftbasind(crbs, init_typ,range{sumtype},sumtype,bsnum,modstruct,min_notes);
         else
              [initind, initday, drxnvec,basind, shftnum,asympvl]=getshftbasind(crbs, init_typ,range{sumtype},sumtype,bsnum,modstruct,min_notes);
         end
          [redinitind, indvls,drxnvecnew,asympvlnew]=getredinitind2(crbs, initind,drxnvec,asympvl,minshift(sumtype));
          initindnew=redinitind;
          initdaynew=initday(indvls);
          shftnumnew=shftnum(indvls);
          datst{1}.ind=initindnew;
          datst{2}.ind=initindnew;
          datst{3}.ind=basind;
          
          datst{1}.day=initdaynew
          datst{1}.asympvl=asympvlnew;
          datst{2}.day=initdaynew
          datst{2}.asympvl=asympvlnew;
          datst{1}.shftnum=shftnumnew;
          datst{2}.shftnum=shftnumnew;
          datst{1}.drxnvec=drxnvecnew;
          datst{2}.drxnvec=drxnvecnew;
     
              
              for runvl=1:3
             	if runvl==2&length(crbs.allnote)>1
                    nt=crbs.allnote(2)
                else
                    nt=crbs.allnote(1);
                end
                %looping through the four categories to fill
                 [datout]= setvals(datout,nt, runvl, crbs,datst,bsnum)
                
              end
         
    end
      
         [sumplot]=mksumplot2(datout);
         [sumplot]=calcoff(sumplot);
         [sumplot]=calcave(sumplot,sumtype);
         if sumtype==1
            sumplotout=sumplot;
            
         else
            sumplotrevout=sumplot;
         end
         
         end


    function [sumplot]=mksumplot2(datout);
            for runvl=1:3
                sumplot(runvl).ac=datout(runvl).acz
                sumplot(runvl).mu=datout(runvl).muz;
                sumplot(runvl).acpre=datout(runvl).acprez;
                sumplot(runvl).acpst=datout(runvl).acpstz;
                sumplot(runvl).acerrz=datout(runvl).acerrz;
                sumplot(runvl).muerrz=datout(runvl).muerrz;
%                 sumplot(runvl).targeff=datout(runvl).targeff;
%                 sumplot(runvl).contreff=datout(runvl).contreff;
                sumplot(runvl).indnum=datout(runvl).ind;
                sumplot(runvl).bsnum=datout(runvl).bsnum;
                sumplot(runvl).shftnum=datout(runvl).shftnum;
                
                if(runvl<3)
                    sumplot(runvl).drxnvec=datout(runvl).drxnvec;
                    sumplot(runvl).initdays=datout(runvl).day;
                    sumplot(runvl).asympvl=datout(runvl).asympvl;
                end
            end
    
        
function [initindcomb,initdaycomb, drxnvec_comb,basind, shftnumcomb,asympcomb]=getshftbasind(crbs,init_typ,range,sumtype,bsnum,modstruct,min_notes);
   trgnt=crbs.allnote(1);
   asympcomb=[];
   initindcomb=[];
   initdaycomb=[];
   shftnumcomb=[];
   drxnvec_comb=[];
   %need to write code for the other option
   if (init_typ==1)
       if sumtype==1
           inruns=crbs.shiftruns
           intimes=crbs.adjshifttms
           indrxn=crbs.drxn
           
          %this code is for over_ride
          if(~isempty('modstruct'))
              if(modstruct.bsnum==bsnum)
                  inruns{modstruct.shftnum}=modstruct.runs;
                  intimes{modstruct.shftnum}=crbs.flrtmvec(modstruct.runs);
                  
              end
          end
       else
           inruns=crbs.revruns
           intimes=crbs.revshifttms
           indrxn=crbs.drxnrev
       end
       
       
       for jj=1:length(inruns)
           if(~isempty(inruns{jj}))
                shftdys=intimes{jj}
                shftrns=inruns{jj}
                
                if(exist('modstruct'))
                    if(modstruct.bsnum~=bsnum)
                        ind=find(ismember(shftdys,range)==1)
                        
                    else
                        ind=1:length(shftrns)
                    end
                else
                    ind=1:length(shftrns)

                end
                initind=shftrns(ind);
%                 modind=find(crbs.acn(trgnt,initind)>min_notes&crbs.mun(trgnt,initind)>min_notes)
%                 initindmod=initind(modind);
                 initindmod=initind;
                initday=shftdys(ind);
                initdaymod=initday(modind);
                
                shftnum=jj*ones(length(modind),1);
                asympvl=crbs.asympvl{jj}*ones(length(modind),1);
                                
                if(indrxn{jj}=='up')
                        drxnvec=ones(length(initind),1);
                    else
                    drxnvec=zeros(length(initind),1);
                end
                drxnvec_comb=[drxnvec_comb drxnvec']
       
                asympcomb=[asympcomb asympvl'];
                initindcomb=[initindcomb initindmod'];
                initdaycomb=[initdaycomb initdaymod'];
                shftnumcomb=[shftnumcomb shftnum'];
             end
        
       end
       basind=crbs.basruns;
     
   end
        
   
function [sumplot]=calcoff(sumplot);
    indup=find(sumplot(1).drxnvec==1)
    inddn=find(sumplot(1).drxnvec==0)
    
    off=sumplot(1).mu-sumplot(1).ac
    sumplot(1).off(indup)=off(indup);
    sumplot(1).pct(indup)=-off(indup)./sumplot(1).ac(indup);
    sumplot(1).off(inddn)=off(inddn);
    sumplot(1).up(indup)=ones(length(indup),1)
    sumplot(1).pct(inddn)=-off(inddn)./sumplot(1).ac(inddn);
    sumplot(1).drxn(inddn)=zeros(length(inddn),1)
    
    indup=find(sumplot(2).ac>0)
    inddn=find(sumplot(2).ac<0)
    
    sumplot(2).off(indup)=sumplot(2).mu(indup)-sumplot(2).ac(indup);
    sumplot(2).off(inddn)=sumplot(2).mu(inddn)-sumplot(2).ac(inddn);
    
    indup=find(sumplot(3).ac>0)
    inddn=find(sumplot(3).ac<0)
    
    sumplot(3).off(indup)=sumplot(3).mu(indup)-sumplot(3).ac(indup);
    sumplot(3).off(inddn)=sumplot(3).mu(inddn)-sumplot(3).ac(inddn);

  function [sumplot]=calcave(sumplot,sumtype);
  basplot=sumplot(3);
 
  bslist=unique(sumplot(1).bsnum);
  for ii=1:length(bslist)
      ind=find(basplot.bsnum==bslist(ii));
      mnoff(ii)=mean(basplot.off(ind));
  end
   if(~isempty(bslist))
    sumplot(4).bsnum=bslist;
    sumplot(4).off=mnoff;
   end
    
 function [outind, redinds,drxnvecnew,asympvlnew]=getredinitind2(crbs, initind,drxnvec,asympvl,minshift)
     nt=crbs.allnote(1);
     threshary=crbs.acz(nt,initind);
     muary=crbs.muz(nt,initind);
     threshind=find(abs(threshary)>minshift);
     redinitind=initind(threshind)
    nonaind=find(isnan(muary)==0)
    nonaind=initind(nonaind);
    [redinitind2]=intersect(redinitind,nonaind);
    outvec=ismember(initind,redinitind2);
    redinds=find(outvec>0);
    outind=initind(redinds);
    drxnvecnew=drxnvec(redinds);
    asympvlnew=asympvl(redinds);
    
     %FUNCTION SETVALS
     function [datout]= setvals(datout,nt, runvl, crbs,datst,bsnum)
                datvls=datst{runvl}.ind;
                bsvls=bsnum*ones(length(crbs.acz(nt, datvls)),1);
                
                datout(runvl).acz=[datout(runvl).acz crbs.acz(nt, datvls)];
                datout(runvl).muz=[datout(runvl).muz crbs.muz(nt,datvls)];
                datout(runvl).acerrz=[datout(runvl).acerrz crbs.acerrz(nt,datvls)];
                datout(runvl).muerrz=[datout(runvl).muerrz crbs.muerrz(nt,datvls)];
                datout(runvl).acprez=[datout(runvl).acprez crbs.acprez(nt,datvls)];
                datout(runvl).acpstz=[datout(runvl).acpstz crbs.acpstz(nt,datvls)];
                datout(runvl).bsnum=[datout(runvl).bsnum; bsvls]
               
                datout(runvl).ind=[datout(runvl).ind datst{runvl}.ind];
                %SET TARGEFF, CONTREFF
                if(runvl<3)
                    datout(runvl).drxnvec=[datout(runvl).drxnvec datst{runvl}.drxnvec]  
                    datout(runvl).day=[datout(runvl).day datst{runvl}.day]
                     datout(runvl).shftnum=[datout(runvl).shftnum datst{runvl}.shftnum];
                     datout(runvl).asympvl=[datout(runvl).asympvl datst{runvl}.asympvl]
                end
         function [datout]=initfields();
             for runvl=1:3
                 datout(runvl).acz=[];
                 datout(runvl).muz=[];
                 datout(runvl).acerrz=[];
                 datout(runvl).muerrz=[];
                 datout(runvl).acprez=[];
                 datout(runvl).acpstz=[];
                 datout(runvl).bsnum=[];
                 datout(runvl).shftnum=[];
                 datout(runvl).ind=[];
                 datout(runvl).drxnvec=[];
                 datout(runvl).day=[];
                 datout(runvl).asympvl=[];
             end