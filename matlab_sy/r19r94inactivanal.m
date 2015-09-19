%inactivation script  8.10.07
%analysis of  effects of inactivation
%first thing to do is just to compare the average note to catch note for
%frequency and amplitude...(bracket amplitude)



%to do this for each batch
clear strcmd;
clear bt;
%dirchange=1;

dir{1}='/doyale/twarren/r19r94/misc'
bt{1}='batch.keep.rand'
dir{2}='/doyale/twarren/r19r94/misc'
bt{2}='batch.keep.rand.dir'



dir{3}='/doyale/twarren/r19r94/acsfin'
bt{3}='batch.keep'
dir{4}='/doyale/twarren/r19r94/acsf8907'
bt{4}='batch.keep'
dir{5}='/doyale/twarren/r19r94/mu80907'
bt{5}='batch.keep'
trigcheck=1;
dir{6}='/doyale/twarren/r19r94/mu81007'
bt{6}='batch.keep'

dir{7}='/doyale/twarren/r19r94/ac81007'
bt{7}='batch.keep.undir.rand'

dir{8}='/doyale/twarren/r19r94/ac81007'
bt{8}='batch.keep.dir'

dir{9}='/doyale/twarren/r19r94/mu81207'
bt{9}='batch.keep.rand'


%bt{2}='batch25.catch.keep';
%bt{3}='batch.keep.rand';    
   



for ii=1:length(bt)
    
strcmd=['cd ' dir{ii}]
eval(strcmd)
%[avne,t,f]=get_avn('batch.keep','e',0.2,0.2,'','','obs0'); 

%compare the frequency of the notes
    btcr=bt{ii}
    tbinshft=0.01;
    NFFT=1024;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
    fbins=[1800,2600];
    save BINS_B NFFT fbins tbinshft
% frequency analysis just for 'b'
    load BINS_B
    NT='e';PRENT='';PSTNT='';

    
       % edges=[6000:75:8000];
    
    fv=findwnote4(bt{ii},NT,PRENT,PSTNT,tbinshft,fbins,NFFT,1,'obs0');
    
    vals{ii}=getvals(fv,1,'TRIG');
    stdv(ii)=std(vals{ii}(:,2))
    num(ii)=length(vals{ii}(:,1))

% strcmd=['save ' bt '.mat ' 'fv '    'tbinshft']
% eval(strcmd)

end

plotting
if(trigcheck)
    edges=1800:50:2600    
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