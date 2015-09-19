%syntax_movwindow
function []=mastersynscript(plotvl);
if (plotvl==3)
    subplot(413)
%     pth{1}='/oriole4/bk13bk12/1106wn_revconting'
    pth{1}='/oriole4/bk13bk12/1105_APVon'
    pth{2}='/oriole4/bk13bk12/1105_acsfon'
    pth{3}='/oriole4/bk13bk12/1105wnon'
    cvl{1}='batch.keep.rand'
    cvl{2}='batch05.keep'
    cvl{3}='batchcomb.keep.rand'
    ptnt='a'
%     cvl{3}='batch07.keep.rand'
 [combvls,ntlabel,outpct]=bk13bk12synanal(pth,cvl,ptnt);
end

if (plotvl==4)
    subplot(4,1,4)
    
    pth{1}='/oriole4/bk13bk12/1106wn_revconting'
    pth{2}='/oriole4/bk13bk12/1107_APVon'
    pth{3}='/oriole4/bk13bk12/1107_ACSF'

    cvl{1}='batch07.rand.keep'
    cvl{2}='batch.keep.rand'
    cvl{3}='batch07.keep.rand'
    ptnt='b'
    [combvls,ntlabel,outpct]=bk13bk12synanal(pth,cvl,ptnt);
end


function [combvls,ntlabel,outpct]=bk13bk12synanal(pth, cvl,ptnt)
WS=40;
PWS=25;
%pharm time

%mkvlsfiles
[combvls,ntlabel,st]=combinevls(pth,cvl,ptnt)
[outpct,stdout,vlsb,stind,ptint]=calcpct(combvls,ntlabel,WS,PWS,ptnt);
figure
    plot(vlsb(:,1),vlsb(:,2),'k.')

figure
plot(combvls(stind:end-(stind-1),1),outpct(1:end),'k.');
hold on;
for ii=1:length(st)
    plot([st(ii) st(ii)], [0 1],'r-');
end

plot(vlsb(ptint:end-(ptint-1),1),stdout(1:end)*10,'r.')


function [combvls,ntlabel,st]=combinevls(pth,cvl,ptnt)
    combvlsa=[]
    combvlsb=[];
    for ii=1:length(pth)
        cmd=['cd ' pth{ii}];
        eval(cmd);
        bt1=cvl{ii};
        % bt2='batch31dir'
        %bt2='batch29d7805;
%         tbinshft=.01
tbinshft=.003
        NFFT=512;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
        if(ptnt=='a')
            fbins=[2500 3500];
            tbinshft=.003
            NFT=512
        else
            fbins=[4500 5500]
            tbinshft=.03
            NFFT=512
        end
        NT{1}='a';PRENT='';PSTNT=''
        NT{2}='b';
        fva{ii}=findwnote4(bt1,'a',PRENT,PSTNT,tbinshft,fbins,NFFT,1,'obs0');
        fvb{ii}=findwnote4(bt1,'b',PRENT,PSTNT,tbinshft,fbins,NFFT,1,'obs0');
        valsaprestim{ii}=getvals_sec(fva{ii},1,'TRIG');
        valsbprestim{ii}=getvals_sec(fvb{ii},1,'TRIG');
        
        if (ptnt=='a')
            st(ii)=valsaprestim{ii}(1,1);
        else
            st(ii)=valsbprestim{ii}(1,1);
        end
        
    end
    for ii=1:length(pth)
       combvlsa=[combvlsa;valsaprestim{ii}]
       combvlsb=[combvlsb;valsbprestim{ii}]
    end
    lna=length(combvlsa(:,1));
    lnb=length(combvlsb(:,1));
    combab=[combvlsa;combvlsb];
    
    [sortout,sortind]=sort(combab(:,1));
    combvls=combab(sortind,:);
    inda=find(sortind<=lna);
    indb=find(sortind>lna);
    ntlabel(inda)='a'
    ntlabel(indb)='b';
    
    uniquetms=unique(combab(:,1));
    for ii=1:length(uniquetms)
        crind=find(combab(:,1)==uniquetms(ii));
        ln=length(crind);
        randsortvls=randperm(ln);
%         randcrind=crind(randsortvls);
%         combab(crind,:)=combab(randcrind,:);
%         ntlabel(crind)=ntlabel(randcrind);
    end
    %for each of unique values, randomize the position of these values
    
 
        
 function [outpct,stdout,vlsb,stind,ptint]=calcpct(combvls,ntlabel,WS,PWS,ptnt)
     
     WSint=floor(WS/2);
     PTint=floor(PWS/2);
     stind=WSint+1;
     ptint=PTint+1;
     vec=WSint+1:length(ntlabel)-WSint;
     for ii=1:length(vec)
           vl=vec(ii);
             inda=find(ntlabel(vl-WSint:vl+WSint)=='a');
             indb=find(ntlabel(vl-WSint:vl+WSint)=='b');
             outpct(ii)=length(inda)/(length(inda)+length(indb));
     end
    if(ptnt=='a')
        ptind=find(ntlabel=='a');
    else
        ptind=find(ntlabel=='b');
    end
    vlsb=combvls(ptind,:);
    vec=PTint+1:(length(ptind)-PTint)
    for ii=1:length(vec)
        vl=vec(ii);
        pitchvec=vlsb(vl-PTint:vl+PTint,2)
        stdout(ii)=std(pitchvec)./mean(pitchvec);
    end
        
             
             
        
    

%loop though file, calculating the percentage.

