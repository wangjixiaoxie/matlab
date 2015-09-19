%create a running average of the stimind
%and the catchind
vls=getvals_sec(fvpt,1,'trig');
windowsize=30;
WS=windowsize;
fac=1/7;
analysis_days=1
outlier_factor=2.5;
REMOVE_OUTLIERS=1;

lnct=length(crctind);
numct=ceil(fac*lnct);
randsamp=rand(lnct,1);
[sortout,ind]=sort(randsamp);
randind=ind(1:numct);
sortout=sort(randind);


rand_crctind=crctind(sortout);
% 
flrdays=unique(floor(vls(:,1)));
flrvls=floor(vls(:,1));
unqdays=unique(flrdays);
for ii=1:length(unqdays)
    crday=unqdays(ii);
    crind=find(flrvls==crday);
    edges(ii,:)=[crind(1) crind(end)]
    


end

%now go through and separate out indices into groups by day
%by definining an edges matrix

%loop through the edges matrix and find the
%relevant catch and fbind for each day
for ii=analysis_days
    ctchsubind=find(rand_crctind>=edges(ii,1)&rand_crctind<=edges(ii,2));
    fbsubind=find(crfbind>=edges(ii,1)&(crfbind<=edges(ii,2)));
    %1 is catch, 2 is stim
    outvls{1}.tm=vls(rand_crctind(ctchsubind),1);
    outvls{1}.pt=vls(rand_crctind(ctchsubind),2);
    outvls{2}.tm=vls(crfbind(fbsubind),1);
    outvls{2}.pt=vls(crfbind(fbsubind),2);
   
        for jj=1:2
            medvl=median(outvls{jj}.pt);
            stdvl=std(outvls{jj}.pt);
            if(REMOVE_OUTLIERS)
            removeind=find((outvls{jj}.pt<(medvl-outlier_factor*stdvl)))
            removeind2=find((outvls{jj}.pt>(medvl+outlier_factor*stdvl)));
            origind=1:length(outvls{jj}.pt)
            [cmbremoveind]=[removeind; removeind2]
            finalind=setdiff(origind,cmbremoveind);
            end
            
           dataout{jj,ii}.tm=filtfilt(ones(1,windowsize)/windowsize,1,outvls{jj}.tm(finalind))
           dataout{jj,ii}.pt=filtfilt(ones(1,windowsize)/windowsize,1,outvls{jj}.pt(finalind))
           
           dataout{jj,ii}.tm=dataout{jj,ii}.tm(WS/2+1:(end-WS/2-1));
            dataout{jj,ii}.pt=dataout{jj,ii}.pt(WS/2+1:(end-WS/2-1));
           
%         
%              dataout{jj,ii}.tm=medfilt1(outvls{jj}.tm(finalind),windowsize)
%             dataout{jj,ii}.pt=medfilt1(outvls{jj}.pt(finalind),windowsize)
% %             dataout{jj,ii}.tm=resample(dataout{jj,ii}.tm(WS/2:end-WS/2),1,WS/2,1);
%              dataout{jj,ii}.pt=resample(dataout{jj,ii}.pt(WS/2:end-WS/2),1,WS/2,1);
        
        end
end
% for ii=1:2
%     dataout{ii}=resample(dataout{ii}(windowsize:end),1,9);
%     
% end