function []=plt_ff(fv,tbins,fbins,NN,sym,TRIG);
%plt_ff(fv,tbins,fbins,NN,sym,TRIG);
% NN is which freq band of fbins to use
% TRIG if this exist and is ==0 will onl yplot the ones that 
%   did not trigger
% if TRIG exists and == 1 will only plot the ones that trig
% if it does not exist then it will plot all

if (~exist('TRIG'))
	TRIG=-1;
end

if (length(tbins)==1)
	Nrow=1;Ncol=1;
elseif (length(tbins)==2)
	Nrow=1;Ncol=2;
elseif (length(tbins)==3)
	Nrow=1;Ncol=3;
elseif (length(tbins)==4)
	Nrow=2;Ncol=2;
else
	Nrow=3;Ncol=ceil(length(tbins)/3);
end


vals=[];
for ii =1:length(fv)
	if (TRIG==-1)
		vals=[vals;fv(ii).mxvals(:,2:end)];
	else
		if (fv(ii).TRIG==TRIG)
			vals=[vals;fv(ii).mxvals(:,2:end)];
		end
	end
end

ax=zeros([4,1]);
for ii = 1:4
	ax(ii)=subplot(Nrow,Ncol,ii);hold on;grid on;
	[b,a]=hist(vals(ii:length(tbins):end,NN),[fbins(NN,1):fbins(NN,2)]);
	plot(a,b./sum(b),sym);
end
linkaxes(ax);
return;
