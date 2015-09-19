%for simple, starting version of this script
%just identify prelesion and postlesion pitch shifts, sequence shifts.

%then for each write a baseline, shift, and recovery value
%write to a matrix for plotting.


function []=lesanal(bs,bnum)
for ind=1:length(bnum)
    crbs=bnum(ind);
    crpth=bs(crbs).path;
    cmd=['cd ' crpth]
    eval(cmd)
    load sumdata.mat
    
    
    LESDATE=datenum(avls.LESIONDATE)
    prelesind=[];
    postlesind=[];
    
    %divide up all the individual experiments into prelesion and postlesion
    for ii=1:length(avls.bastms)
        if(datenum(avls.bastms{ii}{1})<LESDATE)
            prelesind=[prelesind ii]
        else
            postlesind=[postlesind ii]
        end
    end
    %now look at the prelesion experiments
    %loop through them and if they are a pitch exp, run ptanal
    %if they are a seq exp, run seqanal.
    
    prelesct=1
    preptct=1;
    pstptct=1;
    pstlesct=1;
    for indct=1:length(prelesind)
        crind=prelesind(indct)
        if(avls.type{crind}=='PT')
            
            outmt(crbs).prelespt(preptct)=ptanal(crind);
            preptct=preptct+1;
            outmt(crbs).drxn=avls.drxn{crind}
        else
            
            outmt(crbs).preles_sq(prelesct)=seqanal(crind);
            prelesct=prelesct+1;
        end
    end
    
    %do the same for the postlesion experimnents
     for indct=1:length(postlesind)
        crind=postlesind(indct)
        if(avls.type{crind}=='PT')
            
            outmt(crbs).pstlespt(pstptct)=ptanal(crind);
            outmt(crbs).drxn=avls.drxn{crind}
            pstptct=pstptct+1;
        else
            
            outmt(crbs).pstles_sq(pstlesct)=seqanal(crind);
            pstlesct=pstlesct+1;
        end
    end
    
 end
   cmd=['save /oriole7/dir1/lesionsumdata/sumdata.mat  outmt']
   eval(cmd);
    
    function [crot]=ptanal(crind)
         load sumdata.mat
         basvls=[];
         crbas=dvl.bas_sl(crind)
         crwn=dvl.wn_sl(crind)
         cr_rec=dvl.rec_sl(crind)
         if(isfield(dvl,'SPLITDAY'))
             if((dvl.SPLITDAY))
                SPLIT=1;
             else
                 SPLIT=0;
             end
         else
                SPLIT=0;
         end
         
         
             for ii=1:length(crbas)
                basvls=[basvls ;vlsout{crbas(ii)}];
             end
             wnvls=vlsout{crwn};

            if(SPLIT)
               pmind=find(mod(basvls(:,1),1)>.5)
               basvls=basvls(pmind,:);
               pmind=find(mod(wnvls(:,1),1)>.5)
               wnvls=wnvls(pmind,:);
            end
             
         crot.freqbas=median(basvls(:,2))
         
         crot.freqsh=median(wnvls(:,2))
%          if (~isempty(cr_rec)&~isnan(cr_rec))
%             crot.freqrec=median(vlsout{cr_rec}(:,2))
%          else
%              crot.freqrec=[NaN];
%          end
         crot.drxn=avls.drxn{crind};
   function [crot]=seqanal(crind)
       load sumdata.mat
       crbas=dvl.bas_sl{crind}
       crwn=dvl.wn_sl(crind);
       cr_rec=dvl.rec_sl(crind);
       
       crot.seqbas=outnotect(crbas,avls.SEQTRGNT)./sum(outnotect(crbas,avls.ALLSEQNT));
       crot.seqsh=outnotect(crwn,avls.SEQTRGNT)./sum(outnotect(crwn,avls.ALLSEQNT));
        if (~isempty(cr_rec)&~isnan(cr_rec))
            crot.seqrec=outnotect(cr_rec,avls.SEQTRGNT)./sum(outnotect(cr_rec,avls.ALLSEQNT))
         else
             crot.seqrec=[NaN];
         end