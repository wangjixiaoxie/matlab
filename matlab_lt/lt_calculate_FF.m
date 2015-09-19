

ADDX=input('ADDX (1 or 0, default 0)? ');
vals = evtaf_freq('batch.catch.keep', [5500 8000], 'b', 128, 'obs0', ADDX);



mean(vals(:,2))
median(vals(:,2))
std(vals(:,2))
