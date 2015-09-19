function [timevals,distnt] = jcrepeatdist(note,batchnotes)
% jcrepeatdist('a')
labels=notestatsJC(batchnotes);
fvals=findwnoteJC(batchnotes,note,'','',0,[2000 2700],8500,1,'obs0',1);
tvals=timing3(fvals);

% labels always starts and ends with '/',
% find places where it goes from not 'a' to 'a'

    isnt=[];
    for i=2:length(labels)
        isnt(i)=isequal(labels(i),note);
    end
    indnt=find(isnt);
    firstnt=find(~isnt(indnt-1));
    lastnt=find(~isnt(indnt+1));
    distnt=lastnt-firstnt+1;
    timevals=tvals(firstnt);