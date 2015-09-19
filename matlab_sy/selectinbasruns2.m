%purpose of this code is to go through sumbs and calculate
%the offz for each initrun
%for each control run
%and for each baseline runs
%save to a masterstruct suminitshs
%with initacz, initmuz, initoff, init_bsvl, 


%this option tells whether to plot initind, selected as consecutive mu runs
%over some threshold.
function [sumplot]=selectinbasruns2(sumbs,init_typ, range, minshift)

%init_typ 1 means days since wn was turned on.
%init_typ 2 means days since crossing threshold

%range is values, for example 1:5, includes first five days since wn was
%turned on.

for jj=1:6
              initary{jj}=[];
              ctrlary{jj}=[];
              vlsary{jj}=[];
              basary{jj}=[];
end
          initdaycomb=[];
          bsinitnumcomb=[];
          initindcomb=[];
          bsbasnumcomb=[];
          basindcomb=[];
          bsctrlnumcomb=[];
          bsctrlindcomb=[];
          shftnumcomb=[];
for ii=1:length(sumbs)
         crbs=sumbs(ii)
         [initind, initday, basind, shftnum]=getshftbasind(crbs, init_typ,range);
         allnote=crbs.allnote;
          acary{1}=crbs.acz
          acary{2}=crbs.muz
          acary{3}=crbs.acerrz
          acary{4}=crbs.muerrz
          acary{5}=crbs.combeff(crbs.ntind,:);
          if(length(allnote)>1)
            acary{6}=crbs.combeff(crbs.allnote(2),:);
          else
              acary{6}=zeros(1,length(crbs.acerrz));
          end
          targnt=allnote(1);
          [redinitind, indvls]=getredinitind(acary, initind,targnt,minshift);
          initindnew=redinitind;
          initdaynew=initday(indvls);
          shftnumnew=shftnum(indvls);
          for ntind=1:length(allnote)
             nt=allnote(ntind);
             
             %looping through the four categories to fill
             for jj=1:6
                 %effvl
                 if(jj<5)
                    vlsary{jj}=acary{jj}(nt,initindnew);
                 else
                     vlsary{jj}=acary{jj}(initindnew)
                 end
                    if(ntind==1)
                        initary{jj}=[initary{jj} vlsary{jj}]
                        
                        if(jj==1)
                            initdaycomb=[initdaycomb initdaynew];  
                            bsnum=ii*ones(length(vlsary{jj}),1);
                            bsinitnumcomb=[bsinitnumcomb bsnum']
                            initindcomb=[initindcomb initindnew];
                            shftnumcomb=[shftnumcomb shftnumnew];
                         end
                     else
                        ctrlary{jj}=[ctrlary{jj} vlsary{jj}]
                        if(jj==1)
                            bsctrlnum=ii*ones(length(vlsary{jj}),1);
                            bsctrlnumcomb=[bsctrlnumcomb bsctrlnum']
                            bsctrlindcomb=[bsctrlindcomb initindnew];
                        end
                    end

             end
         end
         nt=allnote(1);
         for jj=1:4
                    vlsary{jj}=acary{jj}(nt,basind);
                    basary{jj}=[basary{jj} vlsary{jj}]
             if jj==1   
                    bsnum=ii*ones(length(vlsary{jj}),1);
                    bsbasnumcomb=[bsbasnumcomb bsnum'];
                    basindcomb=[basindcomb basind];
             end               
         end
end
             
[sumplot]=mksumplot(initary,basary,ctrlary, initdaycomb, bsinitnumcomb, initindcomb, bsbasnumcomb, basindcomb,bsctrlnumcomb, bsctrlindcomb, shftnumcomb);
[sumplot]=calcoff(sumplot);
[sumplot]=calcave(sumplot);
function [sumplot]=mksumplot(initary,basary,ctrlary,initdaycomb, bsinitnumcomb, initindcomb, bsbasnumcomb, basindcomb, bsctrlnumcomb,bsctrlindcomb,shftnumcomb)

    %shiftdays
    sumplot(1).ac=initary{1}
    sumplot(1).mu=initary{2}
    sumplot(1).acer=initary{3}
    sumplot(1).muer=initary{4}
    sumplot(1).targeff=initary{5}
    sumplot(1).contreff=initary{6};
    %ctrl days
    sumplot(2).ac=ctrlary{1}
    sumplot(2).mu=ctrlary{2}
    sumplot(2).acer=ctrlary{3}
    sumplot(2).muer=ctrlary{4}
    sumplot(2).targeff=ctrlary{5}
    sumplot(2).contreff=ctrlary{6}
    %bas days
    sumplot(3).ac=basary{1}
    sumplot(3).mu=basary{2}
    sumplot(3).acer=basary{3}
    sumplot(3).muer=basary{4};
    sumplot(3).targeff=initary{5}
    sumplot(3).contreff=initary{6}
    
    
    sumplot(1).initdays=initdaycomb;
    sumplot(1).bsnum=bsinitnumcomb;
    sumplot(1).shftnum=shftnumcomb;
    sumplot(2).shftnum=shftnumcomb;
    sumplot(3).shftnum=shftnumcomb;
    sumplot(2).bsnum=bsctrlnumcomb;
    sumplot(1).indnum=initindcomb;
    
    sumplot(2).indnum=bsctrlindcomb;
    sumplot(3).bsnum=bsbasnumcomb;
    sumplot(3).indnum=basindcomb;
        
function [initindcomb,initdaycomb, basind, shftnumcomb]=getshftbasind(crbs,init_typ,range);
   
   initindcomb=[];
   initdaycomb=[];
   shftnumcomb=[];
   %need to write code for the other option
   if (init_typ==1)
        for jj=1:length(crbs.shiftruns)
            shftdys=crbs.adjshifttms{jj}
            shftrns=crbs.shiftruns{jj}
            ind=find(ismember(shftdys,range)==1)
            initind=shftrns(ind);
            initday=shftdys(ind);
            shftnum=jj*ones(length(initind),1);
            if(~isempty(initind))
                initindcomb=[initindcomb initind'];
                initdaycomb=[initdaycomb initday'];
                shftnumcomb=[shftnumcomb shftnum'];
            end
            
        end
        basind=crbs.basruns;
   end
        
   
function [sumplot]=calcoff(sumplot);
    indup=find(sumplot(1).ac>0)
    inddn=find(sumplot(1).ac<0)
    
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

function [sumplot]=calcave(sumplot);
  basplot=sumplot(3);
  bslist=unique(basplot.bsnum);
  for ii=1:length(bslist)
      ind=find(basplot.bsnum==bslist(ii));
      mnoff(ii)=mean(basplot.off(ind));
  end
  sumplot(4).bsnum=bslist;
  sumplot(4).off=mnoff;
    
 function [redinitind, threshind]=getredinitind(acary, initind,nt,minshift)
     vlsary=acary{1}(nt,initind);
     muary=acary{2}(nt,initind);
     threshind=find(abs(vlsary)>minshift);
     redinitind=initind(threshind);
    nonaind=find(isnan(muary)==0)
    nonaind=initind(nonaind);
    [redinitind,redinds]=intersect(redinitind,nonaind);
    threshind=threshind(redinds);