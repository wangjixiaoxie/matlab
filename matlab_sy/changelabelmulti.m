%changelabelmulti.m
%3.24.08  This script is used to change the labels of multiple batch files
%across directories.  uses the avls.pvls, avls.cvl structure.

orignote='a';
newnote='a-b-cc'


for ii=1:length(avls.pvls)
        strcmd=['cd '  avls.pvls{ii}];
        eval(strcmd);
        bt=avls.cvl{ii}
        changelabel(bt,'a','a-b-cc')
end