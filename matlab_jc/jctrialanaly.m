% jctrialanaly

[vals92r,trigs92r]=triglabel('batch1216','a',1,1,0,0);
allescapes92r=getescapes(trigs92r);
[songnum92r,songpos92r,notetime92r]=jctiming; %note files, then rec files

fvals92r=findwnoteJC('batch1216note','a','','',0,[2000 2700],8500,1,'obs0',0);
for i=1:length(fvals92r)
    shifted92r(i,:)=fvals92r(i).datt;
end
[pitch92r]=jc_pitchmat1024(shifted92r,1024,1020,2,2000,2700,3,'obs0',1);

allescapes92r1=[allescapes92r(1:144) allescapes92r(146:1175)];
FFvals1=mean(pitch92r(800:900,:));


% same position in the next song
count=0;
clear Fe
clear Fesc
for i=1:800
    if allescapes50r(i)==0
        if ~isempty(find(songpos50r(i+1:length(songpos50r))==songpos50r(i)))
            count=count+1;
            Fe(count)=FFvalsR(i);
            Fesc(count)=FFvalsR(i+min(find(songpos50r(i+1:length(songpos50r))==songpos50r(i))));
        end
    end
end
%%%% offset pos
toff50=[];
for ii=1:length(trigs50r)
toff50=[toff50;trigs50r(ii).toffset];
end
toffset50=((toff50/1000)*(32000)-512)/4+240;


%%%%%
ind1=find(allescapes92r==1);
ind0=find(allescapes92r==0);
for i=100:700
    inder=max(find(ind1<i));
    FFind1=ind1(inder-5:inder);
    for jj=1:length(FFind1)
        FFind0(jj)=ind0(max(find(ind0<FFind1(jj))));
        if jj>1 && FFind0(jj)==FFind0(jj-1)
            FFind0(jj-1)=ind0(max(find(ind0<FFind1(jj)))-1);
        end
    end
    hitpredict(i)=median(FFvals92R(FFind0));
    escpredict(i)=median(FFvals92R(FFind1));
    actual(i)=median(FFvals92R(i:i+20));
end
    
