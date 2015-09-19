function AddTemplatesToEvConfig(EvConfigFile,varargin);
% AddTemplatesToEvConfig(EvConfigFile,TemplFile1,...,TemplFileN);
%
%  if a config file is missing the template data this will add it in
%  the first input is the name of the config file the second input is the
%  name of the template file for the first note templates, third input is
%  the template file fo the second note templates ....
%

[ND,OP]=ReadEvTAFv4ConfigFile(EvConfigFile,0);

if (length(ND)~=size(varargin,2))
    disp('You need to include 1 template file for each note template');
    disp(['There are ',num2str(length(ND)),' notes being targetting']);
    disp(['but you only listed ',num2str(size(varargin,2)),' template files']);
    disp('Cannot continue');
    return;
end

for ii=1:size(varargin,2)
    TemplFile=varargin{ii};
    templs=load(TemplFile);
    NTempl=size(templs,2);
    ND(ii).Templ=templs;
    [pth,fn] = fileparts(TemplFile);
    ND(ii).TemplFile=fn;

    if (length(ND(ii).CntRng) ~= NTempl)
        disp(['On note ',num2str(ii),' Template size and count range size do not match']);
        disp(['You need to make sure to set the proper values in uievtafsim']);
        if (length(ND(ii).CntRng) > NTempl)
            disp('Reducing CntRng to same size');
            ND(ii).CntRng=ND(ii).CntRng(1:NTempl);
        elseif (length(ND(ii).CntRng) < NTempl)
            disp('Increasing CntRng to fit - will just copy last element');
            for kk=length(ND(ii).CntRng)+1:NTempl
                ND(ii).CntRng(kk)=ND(ii).CntRng(kk-1);
            end
        end
    end
end
disp(['Rewriting config file named : ',EvConfigFile]);
WriteEvTAFv4ConfigFile(EvConfigFile,ND,OP);
return;