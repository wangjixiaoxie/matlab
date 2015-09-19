clear pera_start


pera_start = datenum('2006-10-23 16:00:00', 'yyyy-mm-dd HH:MM:SS')
pera_end=datenum('2006-10-24 16:00:00','yyyy-mm-dd HH:MM:SS');
perb_start=datenum('2006-10-25 13:00:00','yyyy-mm-dd HH:MM:SS');
perb_end=datenum('2006-10-26 13:00:00','yyyy-mm-dd HH:MM:SS');

%trigvec=[];
%for ii=1:length(fvmas3)
 %   trigvec=[trigvec fvmas3(ii).TRIG];
%end


vals=getvals(fvmasc,2,'TRIG');

indtrig=find(vals(:,3)==1);
%indnotrig=find(vals(:,3)==0);



%indnotriga=find(vals(indnotrig,1)>pera_start&vals(indnotrig,1)<pera_end);
%indnotrigb=find(vals(indnotrig,1)>perb_start&vals(indnotrig,1)<perb_end)
indtriga=find(vals(:,1)>pera_start&vals(:,1)<pera_end);
indtrigb=find(vals(:,1)>perb_start&vals(:,1)<perb_end)
%size(indnotriga)
%size(indnotrigb)

edges=[700:10:900]

trigahist=histc(vals(indtriga,2)/2,edges);
%notrigahist=histc(vals(indnotrig(indnotriga),2)/2,edges);
trigbhist=histc(vals(indtrigb,2)/2,edges);
%notrigbhist=histc(vals(indnotrig(indnotrigb),2)/2,edges);
figure;
%stairs(edges, notrigahist, 'r--');hold on;
stairs(edges,trigahist,'r');
hold on;
%stairs(edges, notrigbhist,'k--');hold on;
stairs(edges,trigbhist,'k')
figure
