function colorplot(Yin,Xin)
% function colorplot(Yin,Xin)
% will plot each row of matrix Yin against each row of Xin,
% one row at a time so that the color cycles each time
colors=['k','b','g','r','c','m','y'];
figure;hold;
for i=1:length(Yin)
	plot(Yin(i,:),Xin(i,:),[colors(1+rem(i,7)),'.'])
end
