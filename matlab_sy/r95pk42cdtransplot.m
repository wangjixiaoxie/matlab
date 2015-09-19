cde{1}='o62bk75_210507_0909.2626.cbin'
cde{2}='o62bk75_210507_0911.2654.cbin'
cde{3}='o62bk75_210507_1025.3365.cbin'
cde{4}='o62bk75_210507_1107.3615.cbin'
cda=cde
cdcd{1}='o62bk75_210507_1502.5294.cbin'

cdeons=[2.728 3.2591 2.338 2.6713]
cdaons=[3.5915 4.0927 3.1582 3.5074] 
cdcdons=[3.4988]

for ii=1:4
    subplot(5,2,2*ii-1)
    [dat1,fs]=ReadCbinFile(cde{ii});
    [dat2,fs]=ReadCbinFile(cda{ii});
    cdebnd=(cdeons(ii));
    evspect(dat1,fs,[1 10000], [cdebnd-.05 cdebnd+.7]);
    if (ii<4)
        axis off
    end
    box off;
    subplot(5,2,2*ii)
    cdabnd=(cdaons(ii));
    evspect(dat1,fs,[1 10000], [cdabnd-.05 cdabnd+.7]);
    box off;
    if (ii<4)
        axis off
    end
end

subplot(5,2,9)
[dat1,fs]=ReadCbinFile(cdcd{1});
evspect(dat1,fs,[1 10000]);
axis off;
box off;
    
cdewn{1}='r95pk42_220507_0659.33.cbin'
cdewn{2}='r95pk42_230507_1548.713.cbin'
cdewn{3}='r95pk42_260507_0714.446.cbin'
cdewn{4}='r95pk42_260507_0820.1152.cbin'

cdewnons=[4.4719   3.2114 2.9654]

