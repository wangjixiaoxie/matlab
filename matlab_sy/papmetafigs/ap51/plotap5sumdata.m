%plotap5_sumscript

figure
% crxvl=[1 2 3 6 7 8]
% mnvlsmod=[0 3.56 -21.8 187.8 96.3 204.2]
% stervls=[0 19 16 34 31.9 41.7]
 
%alternative without flipping
mnvlsmod=[0 33.96 -21.8 187.8 96.3 204.2]
% stervls=[0 8.68 16 34 31.9 41.7]


cvx=[2 7]
cvred=[36.3 30.7]
cverr=[9.9 10.2]

subplot(4,1,1:3) 

for ii=1:length(crxvl)
                    bar(crxvl,mnvlsmod,0.6,'EdgeColor','k','FaceColor','none');
                    hold on;
                    plot([crxvl; crxvl],[[(mnvlsmod)-stervls];[(mnvlsmod)+stervls]],'k')
%                     text(crxvl,mnvls*1.5,['n=' num2str(lnvls)]);
                 
end

subplot(414)
for ii=1:length(cvx)
    bar(cvx,-cvred,0.2,'EdgeColor','k','FaceColor','none');
    hold on;
    plot([cvx; cvx],[[-cvred-cverr];[-cvred+cverr]],'k')
end