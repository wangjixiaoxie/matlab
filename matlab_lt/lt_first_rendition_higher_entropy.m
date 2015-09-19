day_num=13;

figure; scatter(1:size(data.b,1),data.b,'o');

data_first=data.b(transition_times{day_num});

figure(gcf); hold on; scatter(transition_times{day_num},data_first,'ro')

mean_first=mean(data_first)

data_rest=data.b;
data_rest(transition_times{day_num})=[];

data_second=data.b(transition_times{day_num}+1);
mean_second=mean(data_second)

data_third=data.b(transition_times{day_num}+2);
mean_third=mean(data_third)


mean_rest=mean(data_rest)

[