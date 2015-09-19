


for ii=1:length(transitionTW)
    crtrns=transitionTW(ii);
    numvl=crtrns.MIN;
    pre=crtrns.pstay(1,numvl);
    wn=crtrns.pstay(2,numvl);
%     pst=crtrns.pstay(3,numvl);
hold on;
    plot(pre,wn,'ko');
end    