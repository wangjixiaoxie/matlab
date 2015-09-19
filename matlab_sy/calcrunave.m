
%This function runs for TYPE='syn' or TYPE= 'pit'
% For SYN
% % INPUT: invals, a list of characters, and intimes, length of times
% % and averages the fractional value of charlist{1} as a fraction of all characters
% in charlisttogether over a window of windowsize
% 
% OUTPUT:runave - running fractional mean
% tpts - time for each of the fractional mean pts.
function calcrunave(avls)
   MKFV=avls.MKFV;
   WINDOWSIZE=avls.WINDOWSIZE
   for ii=1:length(MKFV) 
       clear runave tpts;
       crind=MKFV(ii);
        crpt=avls.pvls{crind}
        crbt=avls.cvl{crind}
        if(isfield(avls,'baspath'))
            cmd=['cd ' avls.baspath crpt]
        else
            cmd=['cd ' crpt]
        end
            eval (cmd);
            cmd=['load ' crbt '.mat'];
            eval(cmd);
  
%make windowsize even
if(mod(WINDOWSIZE,2)==1)
    w_size=WINDOWSIZE+1
else
    w_size=WINDOWSIZE
end

   %target character is first character in charlist
   TARG_CHAR=avls.NT{1};
    
    %get all the inds which have relevant characters
    combinds=[]
%     for charind=1:length(avls.NT)
%        crchar=avls.NT{charind}; 
%        inds=find(ismember(outstruct{1},crchar));
%        combinds(inds)=;
%     end
    
    %step through these inds, starting windowsize/2;
    start_time=w_size/2;
    
    ctvl=1;
    for crindct=start_time:w_size:length(outstruct{1})-start_time
       
        crindvl=crindct;
     crinds=outstruct{1}(crindct-start_time+1:crindct+start_time);
       
       crtms=outstruct{2}(crindct-start_time+1:crindct+start_time);
       lncharvals=length(find(ismember(crinds,TARG_CHAR)));
       runave(ctvl)=lncharvals/length(crinds);
       tpts(ctvl)=mean(crtms);
       ctvl=ctvl+1;
     
    end
       cmd=['save -append  ' crbt '.mat runave tpts'];
      eval(cmd);   
        
   end
    
    
        
  
    
