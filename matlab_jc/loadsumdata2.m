function [avlsout,graphvalsout]=loadsumdata2(bs, birdind);
pth=bs(birdind).path;
nm=bs(birdind).matfilename
strcmd=['cd ' pth ]
eval(strcmd)
strcmd=['load ' nm]
eval(strcmd)
avlsout=avls
if (exist('graphvals'))
graphvalsout=graphvals
else
    graphvalsout=[];
end