function [val,ind] = findpeaks_JC(in,thr,upper)
% function [val,ind] = findpeaks(in,thr,upper)

val = [];
ind = [];

% stupid method
len = length(in);
ups = []; % positive derivative threshold crosses go here
downs = []; % negative derivative threshold crosses here

% get initial status
if in(1) >= thr
  status = 1;
else
  status = 0;
end

for s = 2:len
  
  if (in(s) >= thr)
    
    if (status == 0) % pos deriv thres cross detected
      ups = [ups s];
      status = 1;
    else % still above threshold
      continue;
    end
    
  else % below threshold
    
    if (status == 1) % neg deriv cross detected
      downs = [downs s];
      status = 0;
    else
      continue;
    end
  
  end
  
end

% now get actual peaks
for j = 1:length(ups)
  
  temp = in(ups(j):downs(j));
  [v,i] = max(temp);
  
  if (v < upper)
    val = [val v];
    ind = [ind ups(j)+i-1];
  end
    
    
end