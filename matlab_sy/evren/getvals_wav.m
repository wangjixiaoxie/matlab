function vals=getvals_wav(fv,NN);
%vals=getvals_wav(fv,NN);
%

vals=zeros([length(fv),1]);
for ii = 1:length(fv)
	vals(ii)=fv(ii).mxvals(NN+1);
end
return;
