function changelabel(bt,orignote, PRENT, PSTNT, newstr,max_onset);
if(~exist('max_onset'))
    max_onset=1000000;
end

prentln=length(PRENT)
fid=fopen(bt,'r');
files=[];cnt=0;
while (1)
	fn=fgetl(fid);
	if (~ischar(fn))
		break;
	end
	cnt=cnt+1;
	files(cnt).fn=fn;
end
fclose(fid);
newstrlen=length(newstr);
origstr=[PRENT orignote PSTNT];
spcnt=0;
for ii=1:length(files)
	fn=files(ii).fn;
	if (exist([fn,'.not.mat'],'file'))
		load([fn,'.not.mat']);
        ind=findstr(labels,origstr);
        ind=ind;
    for jj=1:length(ind)
         if length(labels)>=ind(jj)+prentln+newstrlen-1;
            if(onsets(ind(jj)+prentln)-onsets(ind(jj)+prentln-1)<max_onset)
                labels(ind(jj)+prentln:ind(jj)+prentln+newstrlen-1)=newstr;
            end
         end
    end
    cmd=['save -append ' fn '.not.mat labels'];
    eval(cmd);
    
    end
end