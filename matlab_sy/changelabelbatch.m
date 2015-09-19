function changelabelbatch(avls,orignote, PRENT, PSTNT, newstr, onsetdiff);
for ii=1:length(avls.pvls)
      if (~isempty(avls.pvls{ii}))
        strcmd=['cd '  avls.pvls{ii}];
        eval(strcmd);
        eval('mkdir savenotmat');
        eval('!cp *not.mat savenotmat')
        if(exist('onsetdiff'))
            changelabel(avls.cvl{ii},orignote, PRENT, PSTNT, newstr, onsetdiff);
        else
            changelabel(avls.cvl{ii},orignote, PRENT, PSTNT, newstr);
        end
      end
end