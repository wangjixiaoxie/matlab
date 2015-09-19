function avls=getavls(bs,sumbsvl)
    pth=bs(sumbsvl).path;
    mtname=bs(sumbsvl).matfilename;
%     cmd=['cd ' pth 'datasum'];
cmd=['cd ' pth ];
eval(cmd);
    cmd=['load ' mtname];
    eval(cmd);
    test=1;