%calcdailymeanstd.

%input is vals.

function [os]=calcdailymeanstd(vals)

floorvls=floor(vals(:,1));
unqfloor=unique(floorvls);
os.mn=[];
os.sd=[];
os.dy=[];
initdy=unqfloor(1);
for ii=1:length(unqfloor)
    ind=find(floorvls==unqfloor(ii));
    os.mn=[os.mn mean(vals(ind,2))];
    os.sd=[os.sd std(vals(ind,2))];
    os.dy=[os.dy unqfloor-initdy]
end