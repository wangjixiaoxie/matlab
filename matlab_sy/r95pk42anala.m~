%master script, load and combine various fvs, from postscreen./50microstim,
%and 10 microstim

% pathvl{1}='/doyale/twarren/r95pk42/templtest/'
% catchvl{1}='batch.catch.rand'
% pathvl{2}='/doyale/twarren/r95pk42/r95pk42screen/'
% catchvl{2}='batch2.keep.rand'
% pathvl{3}='/doyale/twarren/r95pk42/wnon/'
% catchvl{3}='batch.keep.catch.rand'


% pathvl{1}='/doyale/twarren/r95pk42/templtest/'
% catchvl{1}='batch.catch.rand'
% pathvl{2}='/doyale/twarren/r95pk42/r95pk42screen/'
% catchvl{2}='batch2.keep.rand'
pathvl{1}='/doyale/twarren/r95pk42/templatest/';
catchvl{1}='batch.keep.rand'
pathvl{2}='/doyale/twarren/r95pk42/templatest2/'
catchvl{2}='batch.keep.rand'
pathvl{3}='/doyale/twarren/r95pk42/templatest3/'
catchvl{3}='batch.keep.catch'

pathvl{4}='/doya/r95pk42/wnona/'
catchvl{4}='batch28.keep.catch'
pathvl{5}='/doya/r95pk42/wnona/'
catchvl{5}='batch30.keep.catch'
pathvl{6}=pathvl{5}
catchvl{6}='batch01.keep.catch'
pathvl{7}=pathvl{6}
catchvl{7}='batch03.keep.catch'
pathvl{8}='/doya/r95pk42/wnoffa/'
catchvl{8}='batchcomb.rand'


usex(3:4)=0;
usex(5:13)=0;



NT{1}='d';
NT{2}='a';

fvcomb=[];

%makefv determines whether frequency and marker information for notesis calculated in fv loop
%below. Information is stored in a <batch>.mat file.  It's necessary to run
%this loop once for syntax analysis below.
makefv=1;

%this determines whether to go through if statement below which controls
%frequency analysis.
freqanlys=0;

calc_prop_note=0;
%stimstart 1/30 12pm

if (makefv)
    for ii=8:8%length(pathvl)
        fvnote={};
        for jj=1:length(NT)
            
            PRENT='';PSTNT='';
            tbinshft=0.005;
            strcmd=['cd '  pathvl{ii}];
            eval(strcmd);
            fvnam{jj}=['fv' NT{jj}]
            bt=catchvl{ii}
            NFFT=1024;%number of data points to FFTstrcmd=strcat('!cd ' dir{i})
            fbins=[2700,3600];
            if(usex(ii))
                strcmd1=['fv{jj}=findwnote4(bt,NT{jj},PRENT,PSTNT,tbinshft,fbins,NFFT,1,''obs0'',1);']
                strcmd2=['[valstrigs,trigs{jj}]=triglabel(bt,NT{jj},1,1,0,1);']
            else
                strcmd1=['fv{jj}=findwnote4(bt,NT{jj},PRENT,PSTNT,tbinshft,fbins,NFFT,1,''obs0'');']
                strcmd2=['[valstrigs,trigs{jj}]=triglabel(bt,NT{jj},1,1,0,0);']
            end
                
                eval(strcmd1);
                %eval(strcmd2);
            fvnote{jj}=NT{jj};
        end
        
            strcmd=['save  ' bt '.mat ' 'fv fvnote tbinshft fbins trigs']
        
           
       
            eval(strcmd)
            
     end
end
            

%2/2 1355 - switch to make template earlier.
%2/5 1635 - switch to make cutoff higher.

%over what range of data will analysis occur?
clear switchdts
bnds{1}='2007-06-24 07:00:00'
bnds{2}='2007-07-07 07:00:00'
switchdt{1}='2007-06-29 7:00:00' %12:00:00'


% dir=[1 1 1 1 0]

for ii=1:length(switchdt)switchdts(ii)=datenum(switchdt{ii},'yyyy-mm-dd HH:MM:SS');
end

%make all the days,
%and then find the switchdt days, and for all those days,
%remove one row, modify the afternoon and morning of the previous.
matrixvals=make_time_indices('2007-06-24','2007-07-07',7,21)
for ii=1:length(switchdts)
%find the indices, the last index is the one where it fits
    ind=find(switchdts(ii)>matrixvals(:,1))
    endct=ind(end)
    matrixvals=[matrixvals(1:ind(endct-1),:);matrixvals(endct+1:length(matrixvals),:)]
    matrixvals(ind(endct-1),2)=switchdts(ii);  matrixvals(ind(endct),1)=switchdts(ii)
end%now what you need to do is to change the ones in between, and you're set.
%loop through every one of the switch days
%calculate the
%overlap thresholds...group points in with next day, for each of these days
%make a figure just showing the first day, and the first change with hit
%rate.
for ii=1:length(bnds)
    bndsjs(ii)=datenum(bnds{ii},'yyyy-mm-dd HH:MM:SS');
end
    fvcomb{1}=[]
    fvcomb{2}=[]
for ii=1:8%length(pathvl)
    
    strcmd=['load ' pathvl{ii} catchvl{ii} '.mat']
    eval(strcmd);
        for jj=1:length(NT)
            fvtmp=fv{jj};
            fvcomb{jj}=[fvcomb{jj} fvtmp];
            
        end
end
vals=getvals(fvcomb{2},1,'TRIG');

fvcomb{1}=transval(fvcomb{1},'d',1);


%analysis of frequencies.
if (freqanlys)
    [mnvl,stdv,htrt]=freqanal(vals,matrixvals)
    figure
    x1=23.75
    x2=30
    y1=3000
    y2=3400;
    ax=fill([x1 x2 x2 x1],[y1 y1 y2 y2],[1 .806 .817]);
    hold on;
    timevals=mean(matrixvals,2)-vals(1,1);
    x=[switchdts-vals(1,1) ;switchdts-vals(1,1)];
    y=[4200 ;4400];
    plot(x,y,'r','Linewidth',3)
    hold on;
    errorbar(timevals, mnvl, stdv,'k+','Linewidth',1)
    x=[switchdts-vals(1,1) ;switchdts-vals(1,1)];
    y=[3600 ;4400];
    plot(x,y,'r','Linewidth',3)
    hold on;
end

figure
plot(vals(indht,1)-vals(1,1), vals(indht,2),'k.','MarkerSize',7);
hold on;
plot(vals(indms,1)-vals(1,1), vals(indms,2),'r.','MarkerSize',7);
    

translist{1}='e'
translist{2}='a'
translist{3}='i'
translist{4}='c'
synstruct=syntaxanal(fvcomb{1},translist);
%%%
            
translist={}
fvcomb{1}
            
            




days=1:length(synstruct)

for ii=1:length(days)
    dv=days(ii);
    efrac(ii)=synstruct{dv}.e/synstruct{dv}.ntrans;
    afrac(ii)=synstruct{dv}.a/synstruct{dv}.ntrans;
    ofrac(ii)=1-efrac(ii)-afrac(ii);
    day(ii)=synstruct{dv}.day
end

figure
x1=3.5
x2=9.7
y1=0
y2=1;
ax=fill([x1 x2 x2 x1],[y1 y1 y2 y2],[1 .806 .817]);
hold on;
plot(day-day(1), efrac,'k-o','Linewidth',2);hold on;
plot(day-day(1),afrac,'b-o')
plot(day-day(1),ofrac,'r-o')
axis([0 14 0 1]);


% plot(days(dayind), bpct(dayind), 'b','Linewidth',3)
% plot(days(dayind),bsepct(dayind),'r','Linewidth',3)
% plot(days(dayind),conpct(dayind),'k','Linewidth',3)
% 
% attl=[]
% endttl=[]
% conttl=[]
% bttl=[]
% 
% %analysis of syntax for a
% for ii=1:length(days)
%     acnt=0;
%     endcnt=0;
%     concnt=0;
%     bcnt=0;
%     dyvl=days(ii);
%     
%     ind=find(vlot==dyvl)
%     for jj=1:length(ind)
%         indvl=ind(jj)
%         if(fvcomb(indvl).trans=='-a')
%             acnt=acnt+1;
%         elseif (fvcomb(indvl).trans=='nd')
%             endcnt=endcnt+1;
%         elseif (fvcomb(indvl).trans=='ii')
%             endcnt=endcnt+1;
%         elseif (fvcomb(indvl).trans(1)=='b'|fvcomb(indvl).trans(2)=='b')
%             bcnt=bcnt+1;
%         else
%             concnt=concnt+1;
%         end
%     end
%     attl(dyvl)=acnt;
%     endttl(dyvl)=endcnt
%     conttl(dyvl)=concnt
%     bttl(dyvl)=bcnt
% end    
%     for ii=1:dyvl
%     
%         sumvl(ii)=bttl(ii)+endttl(ii)+conttl(ii)+attl(ii)
%         if(sumvl(ii)==0)
%             apct(ii)=0
%             endpct(ii)=0
%             conpct(ii)=0
%             bpct(ii)=0
%         else
%         bpct(ii)=bttl(ii)/sumvl(ii)
%         endpct(ii)=endttl(ii)/sumvl(ii)
%         conpct(ii)=conttl(ii)/sumvl(ii)
%         apct(ii)=attl(ii)/sumvl(ii)
%         end
%     end
% end
% dayind=[3:8 11:27 29 33 36 39 42]
% figure
% 
% fill([10 25 25 10],[0 0 1 1],[1 .806 .817])
% 
% hold on
% plot(days(dayind),apct(dayind),'k','Linewidth',3)
% plot(days(dayind),bpct(dayind),'g','Linewidth',3)
% plot(days(dayind),conpct(dayind),'b','Linewidth',3)
% plot(days(dayind),endpct(dayind),'r','Linewidth',3)
% 
% 
%     
% 
% 
% 
