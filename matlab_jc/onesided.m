function abc=onesided(normshift,ons,offs)

% Takes a normalized pitch shift curve or predicted pitch shift curve and
% turns it into a one-sided decay by averaging on both sides of the peak.


abc=zeros(1,1000);
notelength=offs-ons+1;
[toss,maxpoint]=max(normshift);

left=maxpoint-1;
right=notelength-maxpoint;
%abc=normshift(n,notelength);
if right>left
    for m=2:left
        abc(m)=mean([normshift(maxpoint-m) normshift(maxpoint+m)]);
    end
    for m=left+1:right
        abc(m)=normshift(maxpoint+m);
    end
else
    for m=2:right
        abc(m)=mean([normshift(maxpoint-m) normshift(maxpoint+m)]);
    end
    for m=right+1:left
        abc(m)=normshift(maxpoint-m);
    end
end

