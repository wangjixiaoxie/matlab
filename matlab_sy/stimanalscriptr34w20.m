
muinds=[1]
clear fvpt valsa

%%%%
%%%following variables deal with the fact that data was collected with two
%%%programs
SOUNDSTR='wav'
STIMSTR='rig'
REPSOUNDSTR='und'
REPSTIMSTR='im'
REPINDS=[1]
%%%

for ii=1:length(muinds)

    crind=muinds(ii);
    pathvl=tvls.pvls{crind}
    if (exist('baspath'))
        cmd=['cd ' baspath pathvl]
    else
        cmd=['cd ' pathvl]
    end
    sndstr{crind}=SOUNDSTR;
    stimstr{crind}=STIMSTR;
    if(ismember(crind,REPINDS))
       sndstr{crind}=REPSOUNDSTR;
        stimstr{crind}=REPSTIMSTR;
    end
    
    eval(cmd);
    bt=tvls.cvl{crind};
    cmd=['load ' bt '.mat']
    eval(cmd);
        tbinshft=-.02;
 
    NFFT=512;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
    fbins=[6800 8200];
    save BINS_B NFFT fbins tbinshft
    % frequency analysis just for 'b'
    load BINS_B
    NT='b';PRENT='';PSTNT='';
    fvst{crind}=findwnote9(bt,NT,PRENT,PSTNT,tbinshft,fbins,NFFT,1,'obs0',0,[150 50], sndstr{crind},stimstr{crind});
%     fv=findwnote4(bt,NT,PRENT,PSTNT,tbinshft,fbins,NFFT,1,'obs0',0)
valsa{crind}=getvals2(fvpt,1,'TRIG');
end



clear notind fbind ctind
for ii=1:length(muinds)
run_num=muinds(ii);
   notind{run_num}=[]
    fbind{run_num}=[];
    ctind{run_num}=[];
    fvt=fvst{run_num};
    bt=tvls.cvl{run_num}
   
    for ii=1:length(fvt)
        if(fvt(ii).STIMTRIG)
            if(fvt(ii).STIMCATCH)
                ctind{run_num}=[ctind{run_num} ii]
            else
                fbind{run_num}=[fbind{run_num} ii]
            end
        else
            notind{run_num}=[notind{run_num} ii]
        end
    end
    crfbind=fbind{run_num};
    crnotind=notind{run_num};
    crctind=ctind{run_num};
    pathvl=tvls.pvls{run_num}
    if (exist('baspath'))
        cmd=['cd ' baspath pathvl]
    else
        cmd=['cd ' pathvl]
    end

    eval(cmd);
   
    
    edges=[6300:50:8000]
    ctvls=valsa{run_num}(crctind,2);
    fbvls=valsa{run_num}(crfbind,2);
    ctmean=mean(ctvls);
    fbmean=mean(fbvls);
    stdfb=std(fbvls);
    stdct=std(ctvls);
    
    hstoutctind=histc(ctvls,edges);
    hsctnrm=hstoutctind./length(hstoutctind);
    hstoutfbind=histc(fbvls,edges);
    hsfbnrm=hstoutfbind./length(hstoutfbind);
    crvls=valsa{run_num};
 
    cmd=['save -append ' bt '.mat crfbind crnotind crctind hsctnrm hsfbnrm crvls ctmean fbmean stdfb stdct']  ;
    eval(cmd);
end

%make hists
% for ii=5:5  2:length(muinds)
% 
%     hstoutctind{ii}=histc(valsa{ii}(ctind{ii},2),edges)
%     hsctnrm{ii}=hstoutctind{ii}./length(hstoutctind{ii})
%     hstoutfdind{ii}=histc(valsa{ii}(fbind{ii},2),edges)
%     hsfbnrm{ii}=hstoutfdind{ii}./length(hstoutfdind{ii})
% 
% hs{1}=histc(valsa{1}(:,2),edges)
% save -append
