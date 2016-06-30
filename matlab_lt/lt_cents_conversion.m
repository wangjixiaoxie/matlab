%% from abs shift to cents
baseline=3000;
shift=72;

C=1200*log2((baseline+shift)/baseline)


%% ==== from cents to absolute shift
cents=-5;
baseline=3000;

final=baseline*2^(cents/1200)
