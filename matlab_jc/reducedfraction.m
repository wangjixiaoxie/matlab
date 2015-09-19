% w22pk23
labels=notestats
labels=labels'
MIN=4;
songstarts=find(labels=='/');
totalnts=0;
targnts=0;
for i=1:length(songstarts)-1
    if ~isempty(find(labels(songstarts(i)+1:songstarts(i+1)-1)~='-')) % is the song labeled?
        totalnts=totalnts+songstarts(i+1)-songstarts(i)-1;
        for k=songstarts(i)+1:songstarts(i+1)-MIN % adjust this
            if isequal(labels(k:k+MIN-1),'aaaa')
                targnts=targnts+1;
            end
        end
    end
end
targnts
totalnts
targnts/totalnts

% #4
% pre=0.0884 (180/2037)
% post (final day) = 0.0181 (24/1323)

% #6
% pre= 0.0691 (288/4166)
% post (final day) = 0.0476 (124/2607)

% #7
% pre= 0.1710 (455/2661)
% post (final day) = 0.0265 (82/3096)
 
% #9
% pre= 0.1989 (747/3755)
% post (final day) = 0.1481 (308/2080)
 
% #10
% pre= 0.0940 (185/1968)
% post (final day) = 0.0227 (19/836)

% #12
% pre= 0.0428 (339/7925)
% post (final day) = 0.0034 (14/4171)

% #16
% pre= 0.1535 (524/3413)
% post (final day) = 0.1169 (207/1770)

% #18
% pre= 0.0292 (203/6956)
% post (final day) = 0.0157 (47/3003)

% #19
% pre= 0.0617 (89/1442)
% post (final day) = 0.0189 (9/475)

