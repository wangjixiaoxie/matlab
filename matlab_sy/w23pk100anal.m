%Xstimanalscript
%analysis of acute effects of X stimulation
%first thing to do is just to compare the average note to catch note for
%frequency and amplitude...(bracket amplitude)

%also want to compare total number of songs [bracket]
%proportion of as in all songs
%distribution of number of repeats

%updated on 12.2.2006, 
    %to make more mechanistically sound.
    %I am going to save fv structure to a datafile,
    %which I am then going to have for eternity



%to do this for each batch
clear strcmd;
clear bt;
%dirchange=1;
%dir{1}='doyale4/twarren/bu92bu1/stim2'
%dir{2}='doyale4/twarren/bu92bu1/stim1'
%dir{3}='doyale4/twarren/bu92bu1'

trigcheck=1;

bt='batch';


%bt{2}='batch25.catch.keep';
%bt{3}='batch.keep.rand';    
    [avn,t,f]=get_avn('batch','b',0.2,0.2,'','','obs0');

clear bt tbinshft NFFT fbins NT PRENT PSTNT

%compare the frequency of the notes
bt='batch.keep.rand'
tbinshft(1)=0.005;
tbinshft(2)=0
NFFT(1)=1024
NFFT(2)=512;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
fbins{1}=[3000,4000];
fbins{2}=[4300 5300]
save BINS_B NFFT fbins tbinshft bt
% frequency analysis just for 'b'
load BINS_B
NT{1}='e';PRENT{1}='e';PSTNT{1}='e';
NT{2}='b';PRENT{2}='b';PSTNT{2}='b';
    
       % edges=[6000:75:8000];
 for ii=1:length(tbinshft)   
    fv{ii}=findwnote4(bt,NT{ii},PRENT{ii},PSTNT{ii},tbinshft(ii),fbins{ii},NFFT(ii),1,'obs0');
    
    vals{ii}=getvals(fv{ii},1,'TRIG');

    
    if ii==1
    strcmd=['save ' bt '.mat ' 'fv '    'tbinshft ' 'NFFT ' 'fbins ' 'PRENT ' 'PSTNT ' ]
    eval(strcmd)
    else
    strcmd=['save -append ' bt '.mat ' 'fv '    'tbinshft ' 'NFFT ' 'fbins ' 'PRENT ' 'PSTNT']
    
    eval(strcmd)
    end
 end
plotting
if(trigcheck)
    edges=4300:50:5200    
    indtrig=find(vals(:,3)==1)
        indnotrig=find(vals(:,3)==0);
        trighist=histc(vals(indtrig,2),edges);
        notrighist=histc(vals(indnotrig,2),edges);
        figure  
        stairs(edges,trighist,'r');
        hold on
        stairs(edges, notrighist,'k--');
else
    col={'r' 'k' 'b'}
    figure
    for i=1:length(bt)
        trighist{i}=histc(vals{i}(:,2),edges);
        stairs(edges,trighist{i},col{i});
        hold on;
    end
end