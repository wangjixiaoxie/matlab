
ind=1;
while (ind<size(c,2))
	len=c(2,ind);

	ind1=ind+1;
	ind2=ind1+len-1;

	xv=c(1,ind1:ind2);
	yv=c(2,ind1:ind2);
	hold on; plot([xv,xv(1)],[yv,yv(1)],[clr,'-']);
	ind=ind2+1;
end
	
