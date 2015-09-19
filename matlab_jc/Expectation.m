function [E,L] = Expectation(X,k,W,M,V)
[n,d] = size(X);
a = (2*pi)^(0.5*d);
S = zeros(1,k);
iV = zeros(d,d,k);
for j=1:k,
if V(:,:,j)==zeros(d,d), V(:,:,j)=ones(d,d)*eps; end
S(j) = sqrt(det(V(:,:,j)));
iV(:,:,j) = inv(V(:,:,j));
end
%*************************************************************************
%*************************** MODIFICATION: LOOP REMOVAL *****************
%*************************************************************************
E = zeros(n,k);
% for i=1:n,
% for j=1:k,
% dXM = X(i,:)'-M(:,j);
% pl = exp(-0.5*dXM'*iV(:,:,j)*dXM)/(a*S(j));
% E(i,j) = W(j)*pl;
% end
% E(i,:) = E(i,:)/sum(E(i,:));
% end
chunkSize = 1000;
howManyFullChunks = fix(n/chunkSize);
numberOfRemainingVectors = n - howManyFullChunks*chunkSize;
pdfValueVector = zeros(n,1);
for j = 1:k
for chunkCounter = 1:howManyFullChunks
modificationRange = (chunkCounter-1)*chunkSize+1:chunkCounter*chunkSize;
meanSubtractedChunk = X(modificationRange,:)'-repmat(M(:,j),1,chunkSize);
pdfValueVectorToExp = diag(-0.5*meanSubtractedChunk'*iV(:,:,j)*meanSubtractedChunk);
pdfValueVector(modificationRange) = exp(pdfValueVectorToExp)/(a*S(j));
end
modificationRange = howManyFullChunks*chunkSize+1:n;
meanSubtractedChunk = X(modificationRange,:)'-repmat(M(:,j),1,length(modificationRange));
pdfValueVectorToExp = diag(-0.5*meanSubtractedChunk'*iV(:,:,j)*meanSubtractedChunk);
pdfValueVector(modificationRange) = exp(pdfValueVectorToExp)/(a*S(j));
E(:,j) = W(j)*pdfValueVector';
end
sumInColumnDirection = sum(E,2);
divisor = repmat(sumInColumnDirection,1,k);
L = sum(log(sum(E,2)));
E = E./divisor;
%**************************************************************************
%*************************** MODIFICATION ENDS HERE ***********************
%**************************************************************************

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% End of Expectation %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%