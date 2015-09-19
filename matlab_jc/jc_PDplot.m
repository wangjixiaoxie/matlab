function jc_PDplot(pitch_data,F_low,F_high)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%PLOTTING%%%%%%%%%%%%%%%%%%%%%%%
%Calculate mean to use for the red hashed line
for i=1:length(pitch_data(1).pitches)
    for j=1:length(pitch_data)
        matr(i,j)=pitch_data(j).pitches(i);
    end
    averaged(i)=mean(matr(i,:));
end
mmean=mean(averaged(100:600));
for i=1:length(pitch_data(1).pitches)
    t(i)=mmean;
end



%Plotting - like pplotter


figure; hold on;subplot(141); hold on;xlim([0 800]);ylim([F_low F_high+1800])
for i=1:10
    plotshiftedD=pitch_data(i).pitches+180*(i);
    %plot(plotshiftedD)
    plot(t+180*i,'LineStyle',':','Color',[1 0 0])
    plot(plotshiftedD,'r');
end
subplot(142); hold on; xlim([0 800])
for i=11:20
    plotshiftedD=pitch_data(i).pitches+180*(i);ylim([F_low+1800 F_high+3600])
    %plot(plotshiftedD)
    plot(t+180*i,'LineStyle',':','Color',[1 0 0])
    plot(plotshiftedD,'r');
end