%recovery notes.m
% 
% 4 recovery birds which meet criterion of 3-7 days with at least maintained shift of 2.5 standard deviations 
% of the original mean.  
% % 
% 
% %bk63w43
% %divide by two
% max shift - 6453
% 20.
% day 1 recov - 37   6793
% day 2 recov - 23   6893
% day 3 recov - 26  6995
% 
% initmean - 6970
clear recbs
recbs(1).bas=3485
recbs(1).sh=3227
recbs(1).rec=[3397 3447 3498]

% 
% %pk20r49
% maxshift- 22 - 3543
% day 1 - recov- 3413
% day 2 - reov -3333
% day 3- reov - 3249
% initmean - 3245
recbs(2).bas=3245
recbs(2).sh=3547
recbs(2).mu=3508
recbs(2).std=59.4
recbs(2).rec=[3413 3333 3249]
% 
% % 
% % %pu57w52
% % max shift -7974
% % day 1- 7667
% % day 2-7330
% % day 3- 7247
% 
% initmean - 7185

recbs(3).bas=3593
recbs(3).sh=3987
recbs(3).mu=4042;
recbs(3).std=47.13
recbs(3).rec=[3834 3665 3624]

% %Evren's two birds - lesion 2.
% mx shift  - 2239
% 
% 
% initmean - 2361
% rec- 132.3
% rec2 - 95.4
% rec3- 73.5
% 
recbs(4).bas=2361
recbs(4).sh=[2361+174]
recbs(4).rec=[2361+132 2361+95 2361+73.5]

% o11bk33

recbs(5).bas =2251
recbs(5).sh=2167
recbs(5).rec=[2192 2214 2224]

% recbs=bs;

% day 2- 2214
% day 3- 2224
% 
% % These birds 
% % are : my two recovery birds pk20r49, and bs(10)- pu56w26,
% % 
% %also check out bk63w43 - natural recovery.
% and two of evren's birds: these are lesion bird-2
% % and o11bk33- second shift.
% % 
% % %first bird
% % 2224 - 4th day recovery.
% % 2141- max shift.
% % initial bas- 2244
% % sd=32
% % 80.6 percent recovery
% % 
% % %next bird 
% % lesion bird -2 , load /oriole5/evrenlesiondata/sumdata
% % 
% % shift =173.8
% %rec= 77.0
% %55.7 percent
% 
% %bs(10)
% recovery 92.0 percent.
% 
% bs(1)
% recovery 100 percent
% 
% 
% so, mean recovery is mean (100 +92+55.7 +80.6)
% 82.1 +/-9.7 %
% 
% 
% 
% 
% 
