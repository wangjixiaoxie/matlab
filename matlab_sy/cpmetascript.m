cmd1=['cd /doyale3/twarren/dipper/pk72r70/wnonpostlesion']
destdir='/doya3/pk72r70/wnonpostlesion'
cmd2=['cd ' destdir]
for ii=[15 16 17 19]
    eval(cmd1)
    cpscript(destdir,70,ii)
    eval(cmd2)
    bt=['batch' num2str(ii)];
    cmd3=['ls *70_' num2str(ii) '*cbin>' bt]
    eval(cmd3)
    cleandir4(bt,1e4,500,6,10,'obs0')
    btdcrd=['batch' num2str(ii) '.dcrd'];
    mk_rmdata(btdcrd,1)
    
end
