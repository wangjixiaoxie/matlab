
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

function [shft_s] = shiftanal2(bs, brdind)
minstd=1.5
outct=1;
for ct=1:length(brdind)
    brdindvl=brdind(ct)
     [avls,graphvals]=loadsumdata(bs,brdindvl);
     
     wn=avls.wn
     ntind=bs(brdindvl).ntind;
     contrind=bs(brdindvl).contrind;
     allind1=[ntind contrind]     
     contrvl=[zeros(length(ntind),1); ones(length(contrind),1)]
     
     clear subind runs diff
     outruns=[];
     for zz=1:length(allind1)
         ii=allind1(zz);
        %each row (jj) is a different shift
        %this loop sets the runs{1}(jj,:)-lim runs, 
        %and runs{2}(jj,:)-all runs
        %also sets the mxshift(jj), drxn{jj}
        for jj=1:length(wn(ii).on(:,1))
                flr_wnonall(jj,:)=floor(wn(ii).on(jj,:))
                flr_wnoffall(jj,:)=floor(wn(ii).off(jj,:));
                indon=find(flr_wnonall(jj,:));
                indoff=find(flr_wnoffall(jj,:))
                flr_wnon{jj}=flr_wnonall(jj,indon);
                flr_wnoff{jj}=flr_wnoffall(jj,indoff);
                
                flrtmvec=floor(avls.adjtimes(:,1));
            %subind{1} are all the inds, subind{2} are the asymp inds.
                subind{1}=find(flrtmvec>=flr_wnon{jj}(end)&flrtmvec<flr_wnoff{jj}(end));    
                subind{2}=find(flrtmvec>=flr_wnon{jj}(1)&flrtmvec<flr_wnoff{jj}(end));    
                subind{3}=1:length(flrtmvec);
                %now find which of these runs are muruns
                for subindvl=1:3
                    crind=subind{subindvl};
                    matchind=find(ismember(crind,avls.mulist));
                    runs=crind(matchind);
                    [sortruns,sortind]=sort(flrtmvec(runs));
                    %only reset outruns if run is a target run
                    if(sortind&contrvl(zz)==0)
                        outruns{subindvl}{jj}=runs(sortind);
                        outrnstmp=runs(sortind);
                    
                        if(subindvl==1)
                            aczvls=avls.acz(ii,outrnstmp);
                            absvls=abs(aczvls);
                            indabs=find(absvls>minstd);
                            outruns{subindvl}{jj}=outrnstmp(indabs);
                        end
                     end

                end
                aczall=avls.acz(ii,outruns{2}{jj});
                [tmp,mxind]=max(abs(aczall));
                mxshift(jj)=aczall(mxind);
                diff(jj,:)=mxshift(jj)-avls.acz(ii,:);
                if(mxshift(jj)>0)
                    drxn{jj}='up'
                else
                    drxn{jj}='do'
                end
      end
               
            %for these runs, 
            %calculate 1.  the tmoffset for each day.
               %2.  the acz, mu mnz, and stdz for each day.
               %3.  the diff in z score.
               %4 the %recovery for each day based on mean.
  %to calculate offset
                
                outnm=outct+zz-1;    
                
                shft_s(outnm).flrtmvec=flrtmvec;
                shft_s(outnm).flr_wnon=flr_wnon;
  %to calculate this/need to switch into mu coordinates.
                shft_s(outnm).flr_wnoff=flr_wnoff;
                shft_s(outnm).acz=avls.acz(ii,:);
                shft_s(outnm).muz=avls.muz(ii,:);
                shft_s(outnm).acerrz=avls.acerracz(ii,:)
                shft_s(outnm).mxshift=mxshift
                shft_s(outnm).muerrz=avls.mustderrz(ii,:)
%                 shft_s(outct).diffz=shft_s(outct).acz-shft_s(outct).muz;
%                 shft_s(outct).pctrcv=shft_s(outct).diffz./shft_s(outct).acz*100;
%                 shft_s(indct).pcterrz=shft_s(indct).differrz./shft_s(indct).acerrz;
                shft_s(outnm).shiftruns=outruns{2}
                shft_s(outnm).subruns=outruns{1};
                shft_s(outnm).allruns=outruns{3}(1,:);
                shft_s(outnm).allnote=allind1;
                
                %if control note
                if(contrvl(zz))
                    shft_s(outnm).ntype='ctrl'
                %if target note
                else
                    shft_s(outnm).ntype='targ'
                    %want to identify the matching control structure ind.
                        [shft_s(outnm).effvls,shft_s(outnm).ac_cv, shft_s(outnm).mu_cv]=get_effvls(avls,allind1);
                        numcontrind=length(allind1)-1;
                        if(numcontrind)
                            shft_s(outnm).contrshsind=outnm+1:1:outnm+numcontrind
                            shft_s(outnm).contrntind=[allind1(2:end)];    
                        else
                        shft_s(outnm).contrind=0;
                        end
                end
                shft_s(outnm).bname=avls.bname;
                
                
                
     end
        outct=outct+length(allind1);
     end
    
