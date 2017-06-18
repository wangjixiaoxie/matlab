function [bigmatrix, presize] = bigmatrix_helper(switches,thisix,propbj,smo)
%makes a big matrix of all renditions aligned by pitch (block x rend)

presize = max(cellfun(@length, switches.preix(thisix)));
postsize = max(cellfun(@length, switches.postix(thisix)));
bigmatrix = nan(length(thisix),presize+postsize); % blocks x songrends

switchtimes = switches.thistranstime(thisix)';

deletelines = zeros(size(bigmatrix,1),1);
for i = 1:length(thisix) % each block
    if ~isnan(switchtimes(i))
        
        thisfill = propbj(switches.preix{thisix(i)});
        bigmatrix(i, presize-length(thisfill)+1:presize) = thisfill;
        thisfill = propbj(switches.postix{thisix(i)});
        bigmatrix(i, presize+1:presize+length(thisfill)) = thisfill;
        %smooth
        bigmatrix(i,:) = smooth_lv(bigmatrix(i,:),smo);
    else
        deletelines(i) = 1;
    end
    
end
bigmatrix(deletelines>0,:)=[];
% if size(bigmatrix,1)==1
%     bigmatrix = bigmatrix';
% end