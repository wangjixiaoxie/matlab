%WHAT analysis PRODUCES THE BATCH.MAT???
%first loop through all the data to make a list of times,
%freq values, and whether stim or catch
%
function cnsld_dynmcs2(avls)
fvpt_com=[];
crct_com=[];
crfb_com=[]
clear otmt 
 edges=[0 9/24 13/24 17/24 1]
  if(isfield(avls,'baspath'))
        cmd=['cd ' avls.baspath ]
        eval(cmd);
  end
    load sumdata.mat
 for ii=1:15
 out_bnds{ii}=[5000 5600]
 end
 for ii=16:25
 out_bnds{ii}=[4800 5400]
 end
 for ii=26:length(unqdys)
 out_bnds{ii}=[4800 5400]
 end
 boundlookup(1:15)=1
 
 
%this loop is just to create the huge data matrix

for ii=1:length(avls.analinds)
    crind=ii
     pathvl=avls.pvls{crind}
     btvl=avls.cvl{crind};
    if(isfield(avls,'baspath'))
        cmd=['cd ' avls.baspath pathvl]
        eval(cmd);
        cmd=['load ' btvl '.mat']
        eval(cmd);
  end
    
    
    %offset index values for the length of the overall matrix;
    if isempty(fvpt_com)
        addind=0;
    else
        addind=length(fvpt_com);
    end
    
    crct_com=[crct_com crctind+addind];
   
    crfb_com=[crfb_com crfbind+addind];
    fvpt_com=[fvpt_com fvpt];
    
    
    
end    


%this loop is to pick out the unique days

%for each of the unique days

%point of this script is for a given set of batch files
%for each day - take the first n, the last n, and the middle n
%catch trials  catchnum

%and then take all the stim trials within the same time interval.

%for each - remove outliers.
%rewrite this to do a linear fit of the data
% catchnum=100;

vls=getvals_sec(fvpt_com,1,'trig');
vlstmfloor=floor(vls(:,1));
unqdys=unique(vlstmfloor);

for ii=1:length(unqdys)
    %these are all the inds from the day
    crind=find(ismember(vlstmfloor,unqdys(ii)));

    %which are catchinds?
    ct_tmpind=intersect(crind,crct_com);
    
    ct_freqvls=vls(ct_tmpind,2);
    
   %which are stiminds??
   [st_tmpind]=intersect(crind,crfb_com);
   
   st_freqvls=vls(st_tmpind,2);
   
   %remove outliers
   selstind=find(st_freqvls>out_bnds{ii}(1)&st_freqvls<(out_bnds{ii}(2)))
   st_tmvls=vls(st_tmpind(selstind),1);
    st_freqvls=vls(st_tmpind(selstind),2);
    
   selctind=find(ct_freqvls>out_bnds{ii}(1)&ct_freqvls<(out_bnds{ii}(2)))
   ct_tmvls=vls(ct_tmpind(selctind),1);
    ct_freqvls=vls(ct_tmpind(selctind),2);
    
   [Pst]=polyfit(mod(st_tmvls,1),st_freqvls,1);
   [Pct]=polyfit(mod(ct_tmvls,1),ct_freqvls,1);
   
   otmt(ii).stslope=Pst(1);
   otmt(ii).ctslope=Pct(1);
   otmt(ii).stint=Pst(2);
   otmt(ii).ctind=Pct(2);
   otmt(ii).stvls=vls(st_tmpind(selstind),1:2);
   otmt(ii).ctvls=vls(ct_tmpind(selctind),1:2);
   
   %now how to downsample points???
   %divide data for each day into NBINS, take median of each
   %freqvl and timevl for each bin...
   clear vlsin
   vlsin{1}(:,1)=st_tmvls
   vlsin{2}(:,1)=ct_tmvls
    vlsin{1}(:,2)=st_freqvls
   vlsin{2}(:,2)=ct_freqvls
   clear vlsout
   
   for vlnum=1:2
       crvlsin=vlsin{vlnum}
       crlen=length(vlsin{vlnum}(:,1));
      
       for edgeind=1:(length(edges)-1)
          crmodvls=mod(crvlsin(:,1),1);
          inds=find(crmodvls>=edges(edgeind)&crmodvls<=edges(edgeind+1));
          
          if(~isempty(inds))
              vlsout{vlnum}{edgeind}=crvlsin(inds,2)
          else
              vlsout{vlnum}{edgeind}=NaN
          end
       end
       
       
   end
   
   otmt(ii).st_vls=vlsout{1}
   otmt(ii).ct_vls=vlsout{2};
   
%    
%     ind_fb,1));
% %    otmt(ii).fb.tmed(2)=median(vls(midind_fb,1));
%    otmt(ii).fb.tmed(3)=median(vls(endind_fb,1));
%    
%    otmt(ii).ct.init.fvls=vls(initind_ct,2);
%    otmt(ii).ct.mid.fvls=vls(midind_ct,2);
%    otmt(ii).ct.end.fvls=vls(endind_ct,2);
%    
%    otmt(ii).ct.fmed(1)=median(vls(initind_ct,2));
%    otmt(ii).ct.fmed(2)=median(vls(midind_ct,2));
%    otmt(ii).ct.fmed(3)=median(vls(endind_ct,2));
%    
%     otmt(ii).ct.tmed(1)=median(vls(initind_ct,1));
%    otmt(ii).ct.tmed(2)=median(vls(midind_ct,1));
%    otmt(ii).ct.tmed(3)=median(vls(endind_ct,1));
%     else
%         otmt(ii)=otmt(ii-1);
%     end
   
   if(isfield(avls,'baspath'))
        cmd=['cd ' avls.baspath ]
        eval(cmd);
   end
end
if(~exist('sumdata.mat'))
cmd=['save sumdata.mat unqdys otmt crct_com crfb_com vls' ]
   eval(cmd);  

else
cmd=['save -append sumdata.mat unqdys otmt crct_com crfb_com vls' ]
   eval(cmd);  
end
   
   
    


%plotting in a separate program...