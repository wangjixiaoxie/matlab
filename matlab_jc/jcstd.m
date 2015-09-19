function outvec=jcstd(invec)
% Guesses the standard deviation for data with outliers:
% For normally distributed data, mean+1sigma=84th prctile and mean-1sigma=16th prctile
% Thus we guess s.d.=(84th prctile-16th prctile)/2
% But this would be pretty messy for small data sets.  So instead use z=+/-0.5sigma
% Thus we guess s.d. =(69th prctile-31st prctile)


outvec=(prctile(invec,69)-prctile(invec,31));
