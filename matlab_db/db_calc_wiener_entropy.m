function [wiener_entropy] = db_calc_wiener_entropy( spec )
%db_calculate_wiener_entropy Calculates winer entropy for a choosen time window
%   Detailed explanation goes here


%creates a subfield, 'entropy', that is a (trial,time) matrix of
%wiener entropy
wiener_entropy = zeros([size(spec,3) size(spec,2)]);

for j = 1:size(spec,3)
    %wiener entropy = geomean(power spectrum)./mean(power spectrum)
    wiener_entropy(j,:) = log(geomean(spec(:,:,j))./mean(spec(:,:,j)));
end



   

end

