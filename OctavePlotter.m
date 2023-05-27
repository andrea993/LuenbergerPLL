system('LPLL_test.exe > output.txt');

out = dlmread('output.txt', ' ');

L = size(out, 1);

dt = 1e-4;
t = (0:L-1)*dt;

figure(1);
plot(t, out');
legend('raw', 'filtered', 'amplitude', 'phase');

