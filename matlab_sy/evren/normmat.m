function out=normmat(in);
% out=normmat(in);
% normalized each column of in to lie between 0 and 1
%

out = 0.0*in;

for ii = 1:size(in,2)
	tmp = in(:,ii);
	tmp = tmp - min(tmp);
	tmp = tmp./max(tmp);
	out(:,ii) = tmp;
end
return;
