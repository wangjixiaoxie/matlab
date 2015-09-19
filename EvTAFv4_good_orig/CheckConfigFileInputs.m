function ERROR=CheckConfigFileInputs(ND,OP);
%
%

NDtemp=ND(ijk);
if (~isfield(NDtemp,'Templ'))
    disp('There are no templates in this configfile');
    disp('Add templates in using the function : '' AddTemplatesToEvConfig '' ');
    disp('Abort');
    ERROR=1;
    return;
end

for ijk=1:length(ND)
    if (size(N

if (~isfield(NDtemp,''))
    for ijk=1:length(ND)
        ND(ijk).
    end
end

if (~isfield(NDtemp,''))
    for ijk=1:length(ND)
        ND(ijk).
    end
end

if (~isfield(NDtemp,''))
    for ijk=1:length(ND)
        ND(ijk).
    end
end

if (~isfield(NDtemp,''))
    for ijk=1:length(ND)
        ND(ijk).
    end
end

if (~isfield(NDtemp,''))
    for ijk=1:length(ND)
        ND(ijk).
    end
end







if (~isfield(OP,'RawSndTH'))
    %put in a value - this is not that important
    OP.RawSndTH=2500;
end

if (~isfield(OP,'FileBufferLeng'))
    disp(['No Predata buffer lenght in the Config File - putting in defulat of 2 seconds']);
    OP.FileBufferLeng=2500;
end

if (~isfield(OP,'BirdName'))
    OP.BirdName='BirdName';
end

if (~isfield(OP,'DIOTriggerPins'))
    OP.DIOTriggerPins='Dev1/port0/line0:7';
end

if (~isfield(OP,'OutputDataDir'))
    OP.OutputDataDir='D:\';
end

if (~isfield(OP,'SilenceT'))
    OP.SilenceT=1.0;
end
if (~isfield(OP,'AIChans'))
    OP.AIChans='Dev1/ai0';
end

if (~isfield(OP,'FS'))
    OP.FS=32;
end

if (~isfield(OP,'MinInVoltage'))
    OP.MinInVoltage=-1.0;
end

if (~isfield(OP,'MaxInVoltage'))
    OP.MaxInVoltage=1.0;
end
if (~isfield(OP,'MinOutVoltage')
    OP.MinInVoltage=-10.0;
end

if (~isfield(OP,'MaxOutVoltage'))
    OP.MaxInVoltage=10.0;
end
if (~isfield(OP,'OutputSoundFileDir'))
    OP.OutputSoundFileDir='D:\';
end

