%lesion example/ bk1bk29; could be used for other birds.
function [datout]=plt_ex(bs)
%first load data from individual bird.
%index for individual bird.
%bs(1) is bk1 bk29
bsnum=1;
%ptprelesion
pt_days{1}={'2010-9-13' '2010-9-17' }
pt_start(1)=1;
%ptpostlesion
pt_days{2}={'2010-9-20' '2010-9-24' }
pt_start(2)=2;
%seqprelsion
seq_days{1}={'2010-9-11' '2010-9-16'}
seq_start(1)=2;
%seqpostlesion
seq_days{2}={'2010-9-23' '2010-9-28'}
seq_start(2)=1;

cmd=['cd ' bs(1).path];
eval(cmd);
load sumdata.mat;


cd /oriole7/dir1/lesionsumdata
load sumdata.mat;

 outdtpre=outmt(bsnum).prelespt
 outdtpst=outmt(bsnum).pstlespt
 sumdat(1).bas=outdtpre.freqbas
 sumdat(1).sh=outdtpre.freqsh
 sumdat(2).bas=outdtpst(end).freqbas
 sumdat(2).sh=outdtpst(end).freqsh(end)

%THIS LOOP PLOTS FREQUENCY
for ii=1:2
    subplot(2,2,ii);
    dyind=getdays(dayout,pt_days{ii});
    zerotime=floor(dayout(dyind(1)+pt_start(ii)));
    for crind=dyind
       xvls=cntpt(crind).tms-zerotime;
       yvls=cntpt(crind).vl;
       ervls=[makerow(cntpt(crind).vl05);makerow(cntpt(crind).vl95)];
       outdtpre=outmt(bsnum).prelespt
     
       if(~isempty(xvls))
        makeplot(xvls,yvls,ervls,sumdat(ii));
       end
       
       hold on;
        
    end
    
    
end
  
%THIS LOOP PLOTS SEQUENCE
       outdtpre=outmt(bsnum).preles_sq
       outdtpst=outmt(bsnum).pstles_sq
       sumdat(1).bas=1-outdtpre.seqbas
       sumdat(1).sh=1-outdtpre.seqsh
       sumdat(2).bas=1-outdtpst.seqbas
       sumdat(2).sh=1-outdtpst.seqsh



for ii=1:2
    subplot(2,2,ii+2);
    dyind=getdays(dayout,seq_days{ii});
    zerotime=floor(dayout(dyind(1)+seq_start(ii)));
    for crind=dyind
       xvls=cntseq(crind).tms-zerotime;
       yvls=1-cntseq(crind).vl;
       ervls=[makerow(1-cntseq(crind).vl05);makerow(1-cntseq(crind).vl95)];
       
       
       if(~isempty(xvls))
        makeplot(xvls,yvls,ervls,sumdat(ii));
       end
       
       hold on;
        
    end
    
    
end



  
    
    %INPUT:
    %daysbnds is 2 strings specifying bnds of days to plot;
    %dayout is a 1-d vector of 
    %RETURNS:daysind, which gives the indices of dayout
    %which are within daysbnds
 function [daysind]=getdays(dayout,daysbnds)
    bnds=datenum(daysbnds);
    [out,daysind]=intersect(dayout,bnds(1):bnds(2));
    
 function [daysind]=makeplot(xvls,yvls,ervls,sumdat)
    TRIXVL=4;
     
     %first make fill function
    xfillvec=[makerow(xvls) makerow(xvls(end:-1:1))]
    yfillvec=[makerow(ervls(1,:)) makerow(ervls(2,end:-1:1))]
    fill(xfillvec,yfillvec,'c','EdgeColor','none');
    hold on;
    plot(xvls,yvls,'k','Linewidth',1);
    
    plot(TRIXVL, sumdat.bas,'k<','MarkerSize',4);
    plot(TRIXVL, sumdat.sh,'r<','MarkerSize',4);
    



