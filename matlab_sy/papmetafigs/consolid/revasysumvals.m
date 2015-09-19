%revasysumvals.m
%this structure is hacked to specifically pull out pharmacology runs where
%there was baseline shift consolid, and then a reversal of learning
%how are you defining shift - first run at the asy level
    %consolid - at least three days at asympote level
%5 of these 6 are without an intermittent manipulation
    
    sv(1).bsnum=1
sv(1).bas=2
sv(1).sh=[19]
sv(1).con=[22]
sv(1).rec=[38 41]
sv(1).fin=[34]
sv(1).NAT=1;

sv(2).bsnum=2;
sv(2).bas=5;
sv(2).sh=[10]
sv(2).con=[21]
sv(2).rec=[24]
sv(2).fin=[32]
sv(2).NAT=0;

sv(3).bsnum=2;
sv(3).bas=32;
sv(3).sh=39;
sv(3).con=43;
sv(3).rec=46;
sv(3).fin=70;  
sv(3).NAT=0;

sv(4).bsnum=3
sv(4).bas=11;
sv(4).sh=13;
sv(4).con=27;
sv(4).rec=29;
sv(4).fin=49;
sv(4).NAT=0;

sv(5).bsnum=3;
sv(5).bas=31;
sv(5).sh=36;
sv(5).con=41;
sv(5).rec=42;
sv(5).fin=49;
sv(5).NAT=0;

sv(6).bsnum=7
sv(6).bas=3;
sv(6).sh=15;
sv(6).con=20;
sv(6).rec=33;
sv(6).fin=[];
sv(6).NAT=0;


sv(7).bsnum=10;
sv(7).bas=11;
sv(7).sh=14;
sv(7).con=20;
sv(7).rec=23;
sv(7).fin=30;
sv(7).NAT=1;