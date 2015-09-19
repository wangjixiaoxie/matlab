%plotap5_sumscript

figure
crxvl=[5 6 7]
crmnvl=[187.8 96.3 204.2]
crer=[34 31.9 41.7]

crxvl=[1 2 3]
crmnvl=[0 34.0 16.5 ]
crer=[0 8.7 8.6]

% p<.05, paired t-test
    for ii=
                    bar(crxvl,mnvlsmod,0.6,'EdgeColor',col,'FaceColor','none');
                    hold on;
                    plot([crxvl crxvl],[(mnvlsmod)-stervls (mnvlsmod)+stervls],'k')
                    text(crxvl,mnvls*1.5,['n=' num2str(lnvls)]);
                 
                 end