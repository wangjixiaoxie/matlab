function [outhist]=plt_ff2(fv,tbins,fbins,sym,NN,df,TRIG,day,HRRNG,SCALEIT);
%[outhist]=plt_ff2(fv,tbinshft,fbins,sym,NN,df,TRIG,day,HRRNG,SCALEIT);
%
% NN is which freq band of fbins to use
% TRIG if this exist and is ==0 will onl yplot the ones that 
%   did not trigger
% if TRIG exists and == 1 will only plot the ones that trig
% if it does not exist then it will plot all
% set TRIG=-1 to plot all as well
% day = pull out only this day of data
% HRRNG = [min,max] use only data from these hours
%SCALEIT=1 (default) scales byt ehharmonic #

if (~exist('df'))
	df=62.5/2;
else
	if (length(df)==0)
		df=62.5/2;
	end
end

if (~exist('NN'))
	NN=3;
else
	if (length(NN)==0)
		NN=3;
	end
end
	

if (~exist('TRIG'))
	TRIG=-1;
else
    if (length(TRIG)<1)
        TRIG=-1;
    end
end

if (~exist('SCALEIT'))    
    SCALEIT=1;
end

if (exist('day'))
    if (length(day)<1)
        clear day;
    end
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

if (exist('day'))
    days=zeros([length(fv),1]);
    for ii = 1:length(fv)
        days(ii)=fv(ii).day;
    end
    pp=find(days==day);
    if (length(pp)<1)
        disp(['no datas for this day']);
        return;
    end
    fv=fv(pp);
end


if (exist('HRRNG'))
    hours=zeros([length(fv),1]);
    for ii = 1:length(fv)
        hours(ii)=fv(ii).hour;
    end
    pp = find((hours>=HRRNG(1))&(hours<=HRRNG(2)));
    if (length(pp)<1)
        disp(['no data for this date and hour rng']);
        return;
    end
    fv=fv(pp);
end

vals=[];
for ii =1:length(fv)
	if (TRIG==-1)
		vals=[vals;fv(ii).mxvals(2:end)];
	else
		if (fv(ii).TRIG==TRIG)
			vals=[vals;fv(ii).mxvals(2:end)];
		end
	end
end
disp(['Number of notes = ',num2str(size(vals,1))]);

ax=zeros([Nrow*Ncol,1]);
for ii = 1:Nrow*Ncol
	ax(ii)=subplot(Nrow,Ncol,ii);hold on;grid on;
	[b,a]=hist(vals(:,NN),[fbins(NN,1):df:fbins(NN,2)]);
    	outhist(ii).a=a;
    	outhist(ii).b=b;
    	if (SCALEIT)
       	 	%plot(a/NN,b./sum(b),sym);
       	 	stairs(a/NN,b./sum(b),sym(1));
    	else	
        	%plot(a,b./sum(b),sym);
        	stairs(a,b./sum(b),sym(1));
    	end
end
%linkaxes(ax);
return;
