
function [W,M,V] = Maximization(X,k,E)
[n,d] = size(X);
V = zeros(d,d,k);
%*************************************************************************
%*************************** MODIFICATION: LOOP REMOVAL *****************
%*************************************************************************
% W = zeros(1,k); M = zeros(d,k);
% for i=1:k, % Compute weights
% for j=1:n,
% W(i) = W(i) + E(j,i);
% M(:,i) = M(:,i) + E(j,i)*X(j,:)';
% end
% M(:,i) = M(:,i)/W(i);
% end
W = sum(E,1);
M = X'*E;
M = M./repmat(W,d,1);
%**************************************************************************
%*************************** MODIFICATION ENDS HERE ***********************
%**************************************************************************
for i=1:k,
for j=1:n,
dXM = X(j,:)'-M(:,i);
V(:,:,i) = V(:,:,i) + E(j,i)*dXM*dXM';
end
V(:,:,i) = V(:,:,i)/W(i);
end
W = W/n;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% End of Maximization %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%