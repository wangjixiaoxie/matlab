% Run this program after notestats and before jc_plotcbin
% 

song_tally=0;
instance=0;
onsets=[0 0];
for i=1:4800
    if o50_labels(i)=='/'
        song_tally=song_tally+1;
    else
        if o50_labels(i)=='q'
            instance=instance+1;
            onsets(instance,1)=song_tally;
            onsets(instance,2)=o50_onsets(i);
        end
    end
end