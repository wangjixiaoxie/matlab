figure

mainpath='/oriole5/bk54w13/'
dirlist={'screen' 'screen' 'canin'  'probin' 'probin'}
fnlist={'bk54w13_140109_0600.1003.cbin' 'bk54w13_140109_0600.1005.cbin' 'bk54w13_180109_0932.3663.cbin'  'bk54w13_270109_0716.-9697.cbin' 'bk54w13_270109_0727.-9569.cbin'}
bndslist=[4.5 8.5; 4 8; 3.7 7.7;3 7; 2.3 6.3]

ln=length(dirlist);

for ii=1:ln
    exsong.ax=subplot(ln, 1, ii)
    exsong.path=[mainpath dirlist{ii}]
    exsong.fn=fnlist{ii}
    exsong.bnds=bndslist(ii,:)
    plotcbin(exsong)
end

