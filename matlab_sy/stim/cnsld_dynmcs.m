
%first loop through all the data to make a list of times,
%freq values, and whether stim or catch
%

fvpt_com=[];
crct_com=[];
crfb_com=[]
clear otmt

%this loop is just to create the huge data matrix
for ii=1:length(avls.analinds)
    crind=avls.analinds(ii);
     pathvl=avls.pvls{crind}
     btvl=avls.cvl{crind};
    if(isfield(avls,'baspath'))
        cmd=['cd ' avls.baspath pathvl]
    else
        cmd=['cd ' pathvl]
    end
    eval(cmd);
    cmd2=['load ' btvl '.mat'];
    eval(cmd2);
    
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

catchnum=100;

vls=getvals_sec(fvpt_com,1,'trig');
vlstmfloor=floor(vls(:,1));
unqdys=unique(floor(vls(:,1)));

for ii=1:length(unqdys)
    %these are all the inds from the day
    crind=find(ismember(vlstmfloor,unqdys(ii)));
%     [sortout,indout]=sort(vls(crind,1));
    %out of these inds which are ctinds?
    %which are catchinds?
    ct_tmpind=intersect(crind,crct_com);
    
    limvls=vls(:,1);
   
    if (length(ct_tmpind)>catchnum)
    initind_ct=ct_tmpind(1:catchnum);
    endind_ct=ct_tmpind(end-catchnum:end);
    midind_ct=ct_tmpind((end/2-catchnum/2):(end/2+catchnum/2));
    
   tmbndsinit=[limvls(initind_ct(1)) limvls(initind_ct(end))]
   tmbndsmid=[limvls(midind_ct(1)) limvls(midind_ct(end))]
   tmbndsend=[limvls(endind_ct(1)) limvls(endind_ct(end))]
    
   %get the stimind in between each of these values
   [fbind]=intersect(crind,crfb_com);
   
   initfb_tmp=find(limvls(fbind)>=tmbndsinit(:,1)&limvls(fbind)<=tmbndsinit(:,2));
   midfb_tmp=find(limvls(fbind)>=tmbndsmid(:,1)&limvls(fbind)<=tmbndsmid(:,2));
   endfb_tmp=find(limvls(fbind)>=tmbndsend(:,1)&limvls(fbind)<=tmbndsend(:,2));
   
   initind_fb=fbind(initfb_tmp);
   midind_fb=fbind(midfb_tmp);
   endind_fb=fbind(endfb_tmp);
   
   otmt(ii).fb.init.fvls=vls(initind_fb,2);
   otmt(ii).fb.mid.fvls=vls(midind_fb,2);
   otmt(ii).fb.end.fvls=vls(endind_fb,2);
   
   otmt(ii).fb.fmed(1)=median(vls(initind_fb,2));
   otmt(ii).fb.fmed(2)=median(vls(midind_fb,2));
   otmt(ii).fb.fmed(3)=median(vls(endind_fb,2));
   
     otmt(ii).fb.tmed(1)=median(vls(initind_fb,1));
   otmt(ii).fb.tmed(2)=median(vls(midind_fb,1));
   otmt(ii).fb.tmed(3)=median(vls(endind_fb,1));
   
   otmt(ii).ct.init.fvls=vls(initind_ct,2);
   otmt(ii).ct.mid.fvls=vls(midind_ct,2);
   otmt(ii).ct.end.fvls=vls(endind_ct,2);
   
   otmt(ii).ct.fmed(1)=median(vls(initind_ct,2));
   otmt(ii).ct.fmed(2)=median(vls(midind_ct,2));
   otmt(ii).ct.fmed(3)=median(vls(endind_ct,2));
   
    otmt(ii).ct.tmed(1)=median(vls(initind_ct,1));
   otmt(ii).ct.tmed(2)=median(vls(midind_ct,1));
   otmt(ii).ct.tmed(3)=median(vls(endind_ct,1));
    else
        otmt(ii)=otmt(ii-1);
    end
   
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