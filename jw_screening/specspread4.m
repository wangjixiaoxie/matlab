function [specmat,matdev]=specspread4(spec_variable,lag)
%[specmat,devmat]=specspread4(spec_variable,lag)
%specspread4 uses the "1 x maxfreq x tensofthousands" variable created by
%specspread3. it will generate the specmat and variability matrix based on whatever
%lag you want.

specmat=zeros(size(spec_variable,2),size(spec_variable,2));
matdev=specmat;
incl=0;
varincl=0;
%varflag=0;
%for debugging..... get rid:
avgest=zeros(1,128);touse=zeros(1,128);scalr=0;

%%NEED TO USE FILEHIST TO PREVENT LOOKUPS INTO WRONG FILES%%
fileend=size(spec_variable,3);
disp([num2str(size(spec_variable,3)),' segments to go.'])
for z=1:fileend
if z/500==fix(z/500);disp([num2str(z),' ',num2str(incl),', ',num2str(varincl)]);end
	for i=1:size(spec_variable,2)
		scalr=spec_variable(1,i,z);
		if z+lag > size(spec_variable,3) | z+lag <= 0
			inclflag=0;
			varflag=0;
			break %goes to next z
		else
			inclflag=1;
			touse=spec_variable(1,:,z+lag);
			if z<.1*size(spec_variable,3)	%excluding 1st 10% from being used in sampling of deviations.
				avgest=touse.*scalr; %so it will drop out of stdev est.
				varflag=0;
			else	
				avgest=specmat(i,:)./incl;
				varflag=1;
			end
			
			specmat(i,:)=specmat(i,:)+(scalr.*touse);
			matdev(i,:)=matdev(i,:)+((scalr.*touse)-avgest).^2;
		end
	end
	if inclflag==1
		incl=incl+1;
	end
	if varflag==1
		varincl=varincl+1;
	end
end

specmat=specmat./incl;
matdev=sqrt(matdev./(varincl-1)); %to use matdev to report stdev, comment out next line. Checked with highflat.. this is a good std estimate.
matdev=matdev./specmat; %this line to use matdev to report CV.
return
