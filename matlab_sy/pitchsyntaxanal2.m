%6.10.08  This preps inactivation data for further analysis and plotting,
%used with inactivanal, requires structure avls, which has many necessary
%inputs
%acmucompar.m, computes inactivation ratios.
function pitchsyntaxanal2(avls);

%this loads in all the data paths, and variables for analysis.
strcmd=['cd ' avls.sumpath 'datasum']
eval(strcmd)
% strcmd=['load ' avls.mtflnm '.mat'];
% eval(strcmd);

if (avls.mkfv)
    for indvl=1:length(avls.mkfv)
       
        ii=avls.mkfv(indvl);
        if(~isempty(avls.pvls{ii}))
            fvnote={};
            for jj=1:length(avls.NT)
                strcmd=['cd '  avls.pvls{ii}];
                eval(strcmd);
                fvnam{jj}=['fv' avls.NT{jj}];
                bt=avls.cvl{ii};  
                fbins=avls.fbins{jj};
                if(avls.changenote)
                    curnote=jj;
                    curpvl=ii;
                    
                    [fbins]=changefreqbins(avls,curnote,curpvl);
                   
                end
                    
%                 if(avls.usex(ii))
%                     strcmd1=['fv{jj}=findwnote4(bt,avls.NT{jj},avls.PRENT{jj},avls.PSTNT{jj},avls.tshft{jj},avls.fbins{jj},avls.NFFT(jj),1,''obs0'',1);'];
% 
%                 else
                    strcmd1=['fv{jj}=findwnote4(bt,avls.NT{jj},avls.PRENT{jj},avls.PSTNT{jj},avls.tshft{jj},fbins,avls.NFFT(jj),1,''obs0'');'];
                %strcmd2=['[valstrigs,trigs{jj}]=triglabel(bt,NT{jj},1,1,0,0);']
%                 end
            
            eval(strcmd1);

%             if(avls.repeatanal(jj))
%                 fv{jj}=repeatanal(fv{jj});
%             end

        fvnote{jj}=avls.NT{jj};
        if(~isempty(avls.exclude))
            if(avls.exclude(jj))
                fv{jj}=exclude_outlrs(fv{jj},2);
            end
         end
            
            end
          strcmd=['save  ' bt '.mat ' 'fv '] ;  
          eval(strcmd);
        end
    end
            
end


% for ii=1:length(avls.bnds)
%     bndsjs(ii)=datenum(avls.bnds{ii},'yyyy-mm-dd HH:MM:SS');
% end



%THIS MAKES FVCOMB WHICH COMBINES ALL DATA INTO A STRUCT
for jj=1:length(avls.NT)
    fvcomb{jj}=[];
    fvcombdir{jj}=[];
end
    for btind=1:length(avls.pvls);
      
        if(~isempty(avls.pvls{btind}))
           
            strcmd=['load ' avls.pvls{btind} avls.cvl{btind} '.mat'];
            eval(strcmd);
            for jj=1:length(avls.NT)
                %do not include directed songs in this master structure.
                ind=find(btind==avls.diron);
                if isempty(ind)
                    fvtmp=fv{jj};
                    btind
                    fvcomb{jj}=[fvcomb{jj} fvtmp];    
                else
                    fvtmp=fv{jj};
                    valstmp=getvals(fv{jj},1,'trig');
                    fvcombdir{jj}=[fvcombdir{jj} fvtmp];
                    avls.dirvls{jj}{btind}=valstmp;
                end
              end
        end
end


%first make a time vector with list of times
[avls.tmvc]=maketimevec(avls);

%now pull out rawtimevals, with no adjustment
[avls.rawtimes]=timeadj(avls.tmvc,0, [1:length(avls.pvls)],0);
%WHAT??
if(~isempty(avls.diron))
    for ii=1:length(avls.diron)
        indvl=avls.diron(ii);
        avls.dirtimes(ii,:)=avls.rawtimes(indvl,:);
        avls.rawtimes(indvl,:)=[0 0];
    end
end

% this pairs up ac and mu trials;
[avls.aclist,avls.mulist]=groupadj(avls.rawtimes,avls.acon, avls.muon)
%stretch out start and end time to 6-22 for daylight savings possibilities.

%this adjust times to take into accoutn dead_time
[avls.rawtimes]=deadtimeadj(avls);

[avls.adjtimes]=timeadj(avls.rawtimes,avls.acoffset, avls.acon,1);
[avls.adjtimes]=timeadj(avls.adjtimes,avls.muoffset, avls.muon,1);
avls.alltimes=avls.adjtimes;
[avls.adjtimes]=timeadj_pretm(avls);
strcmd=['cd ' avls.sumpath 'datasum'];
eval(strcmd);
avls.fvcomb=fvcomb;
avls.fvcombdir=fvcombdir;
[avls]=findbastimes(avls)

strcmd=['save ' avls.mtflnm '-analdata.mat avls'];
    eval(strcmd);

    
function [exfreqbins]=changefreqbins(avls,curnote,curpvl);
    indnt=find(avls.changenote==curnote);
    if(~isempty(indnt))
        indpvl=find(avls.changeruns==curpvl);
        if(~isempty(indpvl))
            exfreqbins=avls.changebins
        else
            exfreqbins=avls.fbins{curnote}
        end
    else
        exfreqbins=avls.fbins{curnote}
    end;
        
    
        
        
        
    