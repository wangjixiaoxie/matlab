function [avlsout,graphvalsout]=loadsumdata(bs, birdind);
pth=bs(birdind).path;
nm=bs(birdind).matfilename
strcmd=['cd ' pth ]
eval(strcmd)
strcmd=['load ' nm]
eval(strcmd)
avlsout=avls
if (exist('graphvals'))
graphvalsout=graphvals
end