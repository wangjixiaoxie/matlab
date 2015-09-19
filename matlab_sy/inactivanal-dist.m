bnd{1}=datenum([graphvals.date{4} ' ' '06:00:00']
bnd{2}=datenum([graphvals.date{4} ' ' graphvals.timon{4}],'yyyy-mm-dd  HH:MM:SS')
bnd{3}=datenum([graphvals.date{4} ' ' graphvals.timoff{4}],'yyyy-mm-dd  HH:MM:SS')
bnd{4}=datenum([graphvals.date{4} ' ' '21:00:00']

offset=1/24;

valsa1=getvals(avls.fvcomb{2},1,'trig')
valsa2=getvals(avls.fvcomb{3},1,'trig')
valsb=getvals(avls.fvcomb{1},1,'trig')

valsbina1{1}=find(valsa1(:,1)