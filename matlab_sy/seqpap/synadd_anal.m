
% %designed to be called after synanalx in order to do additonal
% analysis.

function []=synadd_anal(avls)
cmd=['cd ' avls.baspath ]
eval(cmd);

cmd=['load sumdata.mat']
eval(cmd);

%divide outnotect into relevand time intervals.

for ii=1:length(avls.bastms)
    bastms=datenum(avls.bastms{ii});
    wntms=datenum(avls.wntms{ii});
    rectms=datenum(avls.rectms{ii});
    basind=find(dayout>=bastms(1) &dayout<=bastms(2));
   
    wnind=find(dayout>=wntms(1) &dayout<=wntms(2));
    recind=find(dayout>=rectms(1) &dayout<=rectms(2));



wnmat(ii).vals=outnotect(wnind,:);
wnmat(ii).days=dayout(wnind)-wntms(1)+1;

basmat(ii).vals=outnotect(wnind,:);
basmat(ii).days=dayout(basind)-bastms(1)+1;

recmat(ii).vals=outnotect(wnind,:);
recmat(ii).days=dayout(recind)-rectms(1)+1;

end