function [specmat,matdev]=specspread5(spec_variable,lag)
%[specmat,devmat]=specspread5(spec_variable,lag)
%specspread5 uses the "1 x maxfreq x tensofthousands" variable created by
%specspread3. it will generate the temporal spread matrix and its cv matrix.
% note that this is a weird kind of cv matrix... how temporal struc varies
%across freqs, rather than how things vary across time.
%will use a special cv estimate-- avgest will be figured first, by pooling
%all freqs.
%ATTEMPT #1 = NO WINDOWING OR OVERLAYING.. RAW OUTPUT! USE SMALL SPEC_VAR!!


specmat=zeros(size(spec_variable,3),size(spec_variable,3));
matdev=specmat;
incl=0;
varincl=0;
%varflag=0;
%for debugging..... get rid:
avgest=zeros(1,size(spec_variable,3));touse=zeros(1,size(spec_variable,3));scalr=0;

disp('Figuring average...')
avgest(1,:)=mean(spec_variable,2);

%%NEED TO USE FILEHIST TO PREVENT LOOKUPS INTO WRONG FILES%%
fileend=size(spec_variable,2);
disp([num2str(size(spec_variable,2)),' freq bins to go. Length in segs: ',num2str(size(spec_variable,3))])
for z=3:fileend %don't really need the first couple bins.. song is filtered anyway.
if z/40==fix(z/40);disp([num2str(z),' ',num2str(incl),', ',num2str(varincl)]);end
	for i=1:size(spec_variable,3)
		scalr=spec_variable(1,z,i);
		if z+lag > size(spec_variable,2) | z+lag <= 2
			inclflag=0;
			varflag=0;
			break %goes to next z
		else
			inclflag=1;
			touse(1,:)=spec_variable(1,z+lag,:);
			%the entire if loop below is from specspread4, avg is figured differently here.
		%	if z<.1*size(spec_variable,3)	%excluding 1st 10% from being used in sampling of deviations.
		%		avgest=touse.*scalr; %so it will drop out of stdev est.
		%		varflag=0;
		%	else	
		%		avgest=specmat(i,:)./incl;
		%		varflag=1;
		%	end
			
			specmat(i,:)=specmat(i,:)+(scalr.*touse);
			matdev(i,:)=matdev(i,:)+((scalr.*touse)-avgest).^2;
		end
	end
	if inclflag==1
		incl=incl+1;
	end
	%more outdated specspr4 stuff:
	%if varflag==1
	%	varincl=varincl+1;
	%end
end
 varincl=incl;

specmat=specmat./incl;
matdev=sqrt(matdev./(varincl-1)); %to use matdev to report stdev, comment out next line. Checked with highflat.. this is a good std estimate.
matdev=matdev./specmat; %this line to use matdev to report CV.
return
