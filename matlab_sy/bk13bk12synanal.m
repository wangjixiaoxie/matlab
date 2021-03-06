%syntax_movwindow
function []=mastersynscript(plotvl);
if (plotvl==3)
    subplot(413)
%     pth{1}='/oriole4/bk13bk12/1106wn_revconting'
    pth{1}='/oriole4/bk13bk12/1105_APVon'
    pth{2}='/oriole4/bk13bk12/1105_acsfon'

    cvl{1}='batch.keep.rand'
    cvl{2}='batch05.keep'
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
%pharm times



%mkvlsfiles
[combvls,ntlabel,st]=combinevls(pth,cvl)
[outpct,stdout,vlsb]=calcpct(combvls,ntlabel,WS);
figure
plot(combvls(WS:end,1),outpct(WS:end),'k.');
hold on;
for ii=1:length(st)
    plot([st(ii) st(ii)], [0 1],'r-');
end

plot(vlsb(WS:end,1),stdout(WS:end)/100,'r.')


function [combvls,ntlabel,st]=combinevls(pth,cvl)
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
        fbins=[2500 3500];
        NT{1}='a';PRENT='';PSTNT=''
        NT{2}='b';
        fva{ii}=findwnote4(bt1,'a',PRENT,PSTNT,tbinshft,fbins,NFFT,1,'obs0');
        fvb{ii}=findwnote4(bt1,'b',PRENT,PSTNT,tbinshft,fbins,NFFT,1,'obs0');
        valsaprestim{ii}=getvals_sec(fva{ii},1,'TRIG');
        valsbprestim{ii}=getvals_sec(fvb{ii},1,'TRIG');
        st(ii)=valsaprestim{ii}(1,1);
        
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
        
 function  [outpct,stdout,vlsb]=calcpct(combvls,ntlabel,WS)
     for ii=WS:length(ntlabel)
             inda=find(ntlabel(ii-(WS-1):ii)=='a');
             indb=find(ntlabel(ii-(WS-1):ii)=='b');
             outpct(ii)=length(inda)/(length(inda)+length(indb));
     end
    indb=find(ntlabel=='a');
    vlsb=combvls(indb,:);
    for ii=WS:length(vlsb)
        pitchvec=vlsb(ii-(WS-1):ii,2)
        stdout(ii)=std(pitchvec);
    end
        
             
             
        
    

%loop though file, calculating the percentage.

