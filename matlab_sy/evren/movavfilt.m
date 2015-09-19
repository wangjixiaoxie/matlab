function [fdat]=movavfilt(vals,win,dt);
%[fdat]=movavfilt(vals,win,dt);

t0=floor(min(vals(:,1)));
t1=ceil(max(vals(:,1)));

if (~exist('win'))
	win=2.0/24.0; % two hour signma
else
	if (length(win)==0)
		win=2.0/24.0; % two hour signma
	end
end
		
if (~exist('dt'))
	dt=0.5*win; % step size
else
	if (length(dt)==0)
		dt=0.5*win; % step size
	end
end

t=t0;fdat=[];
while (1)
	if (t>=t1)
		break;
	end

	pp=find(abs(vals(:,1)-t)<=win);
	if (length(pp)==0)
		t=t+dt;
		continue;
	else
		tmp  = mean(vals(pp,2));
		tmp2 = std(vals(pp,2));
		fdat = [fdat;t,tmp,tmp2];
	end

	t=t+dt;
end
return;
