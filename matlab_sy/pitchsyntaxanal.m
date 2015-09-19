function pitchsyntaxanal(avls);

%this loads in all the data paths, and variables for analysis.
strcmd=['cd ' avls.sumpath 'datasum']
eval(strcmd)
% strcmd=['load ' avls.mtflnm '.mat'];
% eval(strcmd);

if (avls.mkfv)
    for indvl=1:length(avls.mkfv)
        ii=avls.mkfv(indvl)
        fvnote={};
        for jj=1:length(avls.NT)
            strcmd=['cd '  avls.pvls{ii}];
            eval(strcmd);
            fvnam{jj}=['fv' avls.NT{jj}]
            bt=avls.cvl{ii}
            
            if(avls.usex(ii))
                strcmd1=['fv{jj}=findwnote4(bt,avls.NT{jj},avls.PRENT{jj},avls.PSTNT{jj},avls.tshft{jj},avls.fbins{jj},avls.NFFT(jj),1,''obs0'',1);']
                %strcmd2=['[valstrigs,trigs{jj}]=triglabel(bt,NT{jj},1,1,0,1);']
            else
                strcmd1=['fv{jj}=findwnote4(bt,avls.NT{jj},avls.PRENT{jj},avls.PSTNT{jj},avls.tshft{jj},avls.fbins{jj},avls.NFFT(jj),1,''obs0'');']
                %strcmd2=['[valstrigs,trigs{jj}]=triglabel(bt,NT{jj},1,1,0,0);']
            end
            
            eval(strcmd1);
%                eval(strcmd2);
            
            if(avls.repeatanal(jj))
                fv{jj}=repeatanal(fv{jj})
            end

fvnote{jj}=avls.NT{jj};
        end
          strcmd=['save  ' bt '.mat ' 'fv ']   
          eval(strcmd)
    end
            
end
for ii=1:length(avls.bnds)
    bndsjs(ii)=datenum(avls.bnds{ii},'yyyy-mm-dd HH:MM:SS');
end
for jj=1:length(avls.NT)
    fvcomb{jj}=[]
end
for ii=1:length(avls.analind)
    btind=avls.analind(ii);
    strcmd=['load ' avls.pvls{btind} avls.cvl{btind} '.mat']
    eval(strcmd);
    for jj=1:length(avls.NT)
        fvtmp=fv{jj};
        fvcomb{jj}=[fvcomb{jj} fvtmp];    
    end
end
strcmd=['cd ' avls.sumpath 'datasum']
eval(strcmd)
avls.fvcomb=fvcomb
strcmd=['save ' avls.mtflnm '-analdata.mat avls'];
    eval(strcmd);


if avls.supanal
    for ii=1:length(avls.bnds)
        switchdts(ii)=datenum(avls.bnds{ii},'yyyy-mm-dd HH:MM:SS');
        matind{ii}=avls.bnds{ii}(1:10)
    end
%     
%     matrixvals=make_time_indices(matind{1},matind{2},7,21)
%     for ii=1:length(switchdts)
% %find the indices, the last index is the one where it fits
%         ind=find(switchdts(ii)>matrixvals(:,1))
%         endct=ind(end)
%         matrixvals=[matrixvals(1:ind(endct-1),:);matrixvals(endct+1:length(matrixvals),:)]
%         matrixvals(ind(endct-1),2)=switchdts(ii);  matrixvals(ind(endct),1)=switchdts(ii)
%     end

    for ii=1:length(avls.pitchind)
        pitchind=avls.pitchind(ii);
        [outvls.mnvl{pitchind},outvls.stdv{pitchind},outvls.htrt{pitchind},outvls.vlot{pitchind}]=freqanal(fvcomb{pitchind})
    end
    for ii=1:length(avls.synind)
        synind=avls.synind(ii);
        fvcomb{1}=transval(fvcomb{1},1);
        outvls.synstruct=syntaxanal(fvcomb{1},avls.translist);
    end


    avls.fvcomb=fvcomb
strcmd=['cd ' avls.sumpath 'datasum']
eval(strcmd)
    strcmd=['save ' avls.datfile '-analdata.mat  outvls avls'];
    eval(strcmd);
end