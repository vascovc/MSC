%% 2.1
clear all
close all
clc

N = 1E6;
x = rand(N,1);
figure(1)
histogram(x)
%+
dx = 0.005;
[M,x_m,bins] = histcounts(x,1/dx);

P = M./(N*dx);
figure(2)
plot(P,'.-')

%% 2.2
clear all
close all
clc

N = 1E6;
x = rand(N,1);

sample_avg = mean(x);
disp(['Sample average - ', num2str(sample_avg)])
disp(['Average theory - ', num2str(1/2)])

sample_var = var(x);
disp(['Sample variance - ', num2str(sample_var)])
disp(['Variance theory- ', num2str(1/12)])


%% 2.3
clear all
close all
clc

N = 1E6;
x = rand(N,1);
y = sqrt(x);
figure(1)
histogram(y);

dx = 0.005;
[M,x_m,bins] = histcounts(x,1/dx);

P = M./(N*dx);
figure(2)
plot(P,'.-')
