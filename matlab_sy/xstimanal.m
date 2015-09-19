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

bt='batch.catch';


%bt{2}='batch25.catch.keep';
%bt{3}='batch.keep.rand';    
    [avnb,t,f]=get_avn('batch26.catch','b',0.2,0.2,'','','obs0');


%compare the frequency of the notes
subplot(211)
bt1='batch_02.keep.catch'
% bt2='batch31dir'
%bt2='batch29d7805;
tbinshft=.034
NFFT=256;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
fbins=[4500 5500];
save BINS_B NFFT fbins tbinshft
% frequency analysis just for 'b'
load BINS_B
NT='b';PRENT='';PSTNT='';

    
       % edges=[6000:75:8000];
    
    fv=findwnote4(bt1,NT,PRENT,PSTNT,tbinshft,fbins,NFFT,1,'obs0');
%     fvd=findwnote4(bt2,NT,PRENT,PSTNT,tbinshft,fbins,NFFT,1,'obs0');
    vlsb2=getvals(fv,1,'TRIG');
    valsd=getvals_sec(fv,1,'TRIG');
strcmd=['save ' bt '.mat ' 'fv '    'tbinshft']
eval(strcmd)

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