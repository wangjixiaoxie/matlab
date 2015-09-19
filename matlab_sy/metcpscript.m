%metcpscript.m
%11-26-07  I'm using this script

smbdir='/doyale3/twarren/dip5/pk80r6/wnoff2/'
smbswitch=['cd ' smbdir]

homedirprewn='/doya2/pk80r6/wnoff2/'
hompreswitch=['cd ' homedirprewn]

homewn='/doya2/pk80r6/wnon2/'
hompostswitch=['cd ' homewn]


for ii=[20:21]

    eval(smbswitch)

    cpscript(homedirprewn,6,ii);
    eval(hompreswitch)
    mkbtcmd=['ls *6_' num2str(ii) '*cbin>batch' num2str(ii)]
    eval(mkbtcmd)
    bt=['batch' num2str(ii)]
    cleandir4(bt,1e4,500,6,10,'obs0')
    mk_rmdata(['batch' num2str(ii) '.dcrd'],1)
    !csh ./rmdata
end

for ii=[25:26]
     eval(smbswitch)

    cpscript(homewn,6,ii);
    eval(hompostswitch)
    mkbtcmd=['ls *6_' num2str(ii) '*cbin>batch' num2str(ii)]
    eval(mkbtcmd)
    bt=['batch' num2str(ii)]
    cleandir4(bt,1e4,500,6,10,'obs0')
    mk_rmdata(['batch' num2str(ii) '.dcrd'],1)
    !csh ./rmdata

end
