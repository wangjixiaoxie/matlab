x = [0 1 0 0];
y = [0 0 1 0];

[cc, lags] = xcov(x, y);

figure; hold on;
plot(lags, cc, 'k');

title('xcov(x, y)');
xlabel('x leads if peak here <--> y leads if peak here');
