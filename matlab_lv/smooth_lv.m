function smovec = smooth_lv(vec,win)

if ~rem(win,2)
    win = win-1;
end

range = (win-1)/2;

smovec = nan(1,length(vec));

for i = range+1:length(vec)-range
    
    smovec(i) = mean(vec(i-range:i+range));
    
end