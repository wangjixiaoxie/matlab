% mu0 = 140;
% sigma0 = 31;
% mu1 = 70;
% n=14;

%% --- not using paired info
Nmax = 100;
PwrAll = [];
for n = 1:Nmax
mu0 = 40;
sigma0 = 31;
mu1 = 20;

pwr = sampsizepwr('t', [mu0, sigma0], mu1, [], n);
PwrAll = [PwrAll pwr];
end

lt_figure; hold on;
plot(1:Nmax, PwrAll, '-');
line([16 16], ylim)

%% --- using paired info
Nmax = 100;
PwrAll = [];
for n = 1:Nmax
    mu0 = 0;
    sigma0 = 31;
    mu1 = 20;
    
    pwr = sampsizepwr('t', [mu0, sigma0], mu1, [], n);
    PwrAll = [PwrAll pwr];
end

lt_figure; hold on;
plot(1:Nmax, PwrAll, '-');
line([16 16], ylim)