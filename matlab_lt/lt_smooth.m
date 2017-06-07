function Ysm = lt_smooth(Y, kernelSD, binsize)

%% default: gaussian kernel
% Y, data to smooth
% kernelSD for gaussian (in same units as binsize)
% binsize (for gaussian)

if ~exist('binsize', 'var');
    binsize = 1; 
end

%% make kernel

           x = -3*kernelSD:binsize:3*kernelSD; % kernel support
            kern = (1/(sqrt(2*pi)*kernelSD)) * exp(-x.^2/(2*kernelSD^2)); % kernel
 
%% 

Ysm = conv(Y, kern, 'same');
            
