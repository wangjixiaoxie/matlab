% using DB's function.

=



population1=1:1000;
population2=50:1050;
func=@(x,y)mean(x)-mean(y);
samplesize=200;
numofboot=20000;
alpha=0.05;

[ CI_median, CI_lo, CI_high, sample_boot ] = lt_db_sample_boot_stat(population1, population2, func, samplesize, numofboot, alpha );

[CI_lo CI_median CI_high]
