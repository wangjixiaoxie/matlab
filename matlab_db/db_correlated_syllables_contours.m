function [ values ] = db_correlated_syllables_contours(multiple_pitch, syllable, day, time_range)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

figure, hold on
title('all pitch contours')

for i = 1:size(multiple_pitch.pc_all.(syllable){day},2)
    if mean(multiple_pitch.pc_all.(syllable){day}(min(time_range):max(time_range),i)) >...
            median(mean(multiple_pitch.pc_all.(syllable){day}(min(time_range):max(time_range),:)))
        plot(multiple_pitch.pc_all.(syllable){day}(:,i), 'b')
        values.high{i} = multiple_pitch.pc_all.(syllable){day}(:,i);
    elseif mean(multiple_pitch.pc_all.(syllable){day}(min(time_range):max(time_range),i)) <...
            median(mean(multiple_pitch.pc_all.(syllable){day}(min(time_range):max(time_range),:)))
        plot(multiple_pitch.pc_all.(syllable){day}(:,i), 'r')
        values.low{i} = multiple_pitch.pc_all.(syllable){day}(:,i);      
    end
end
hold off

values.high = cell2mat(values.high);
values.low = cell2mat(values.low);

figure, hold on
title('mean pitch contours')

plot(mean(values.high,2),'b')
plot(mean(values.low,2),'r')

legend('above median', 'below median','Location','best')

end

