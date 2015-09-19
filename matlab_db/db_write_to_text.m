mkdir('pitch_exp')
for i = 1:length(FF.a.time_and_FF)
    name = ['bk59bk42_day' num2str(i) '.txt'];
    
    fid = fopen(['pitch_exp/' name],'w');
    display(['day_' num2str(i)])
    for j = 1:length(FF.a.time_and_FF{i})
        fprintf(fid, '%.6f\t%.2f\n',FF.a.time_and_FF{i}(j,:));
        display(['day ' num2str(i) '  syllable ' num2str(j)])
    end
end