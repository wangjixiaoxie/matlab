function [FFstatsfront,FFstatsend,toffsetfront,toffsetend]=summary_statsNP(templa,batch,cntrngfront,cntrngend)
fs=32000;
mk_tempf(batch,templa,2,'obs0');

%%%%% FRONT
    get_trigt2(batch,cntrngfront,0.2,128,1,1);
% targeting
    
    [vals,trigs]=triglabel(batch,'a',1,1,0,1);
    toff=[];
    for ii=1:length(trigs)
    toff=[toff;trigs(ii).toffset];
    end
    stdtoff=std(toff);
    toffsetfront1=240+((toff/1000)*(fs)-512)/4;
    [vals,trigs]=triglabel(batch,'b',1,1,0,1);
    toff=[];
    for ii=1:length(trigs)
    toff=[toff;trigs(ii).toffset];
    end
    stdtoff=std(toff);
    toffsetfront2=240+((toff/1000)*(fs)-512)/4;
    toffsetfront=[toffsetfront1;toffsetfront2];
%%% FFstats
    vvals=evtaf_freq('batch',[2000 2700],'a',128,'obs0',1,1);
    FFstatsfront1=vvals(:,2);
    vvals=evtaf_freq('batch',[2000 2700],'b',128,'obs0',1,1);
    FFstatsfront2=vvals(:,2);
    FFstatsfront=[FFstatsfront1;FFstatsfront2];
%%%%% END
get_trigt2(batch,cntrngend,0.2,128,1,1);

    [vals,trigs]=triglabel(batch,'a',1,1,0,1);
    toff=[];
    for ii=1:length(trigs)
    toff=[toff;trigs(ii).toffset];
    end
    stdtoff=std(toff);
    toffsetend1=240+((toff/1000)*(fs)-512)/4;
    [vals,trigs]=triglabel(batch,'b',1,1,0,1);
    toff=[];
    for ii=1:length(trigs)
    toff=[toff;trigs(ii).toffset];
    end
    stdtoff=std(toff);
    toffsetend2=240+((toff/1000)*(fs)-512)/4;
    toffsetend=[toffsetend1;toffsetend2];

vvals=evtaf_freq('batch',[2000 2700],'a',128,'obs0',1,1);
FFstatsend1=vvals(:,2);
vvals=evtaf_freq('batch',[2000 2700],'b',128,'obs0',1,1);
FFstatsend2=vvals(:,2);
FFstatsend=[FFstatsend1;FFstatsend2];

