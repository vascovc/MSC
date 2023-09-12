
%% 1.1
% considerando que todos os numeros apresentam a mesma probabilidade, 1/M
% ao somarem-se as probs \sum_{i=1}^{M}\frac{x}{M} =
% \frac{1}{M}\sum_{x=1}^{M}x = \frac{1}{M}\frac{M(M+1)}{2}

% para a variancia
% \sum_{x=1}^{M}\frac{(x-<x>^2)^2}{M}=\sum_{x=1}^{M}\frac{x^2-<x>^2}{M}
% = \frac{1}{M}(\sum_{x=1}^{M}x^2-\sum_{x=1}^{M}<x>^2)
% 
%% 1.2
clear all
close all
clc

N = [10E2 10E3 10E4 10E6 10E7];
M = 100;
i = 1;

sample_avg = nan([1 length(N)]);
sample_variance = nan([1 length(N)]);
theory_avg = nan([1 length(N)]);
theory_variance = nan([1 length(N)]);
prob = nan([1 length(N)]);

for N_i=N
x = 1+(M-1)*rand(N_i,1);

sample_avg(i) = mean(x);
sample_variance(i) = var(x);

theory_avg(i) = (M+1)/2;
theory_variance(i) = (M^2-1)/12;

prob(i) = length(x(x<60))./N_i;
i = i+1;
end
figure(1)
plot(log10(N),theory_avg)
hold on
plot(log10(N),sample_avg,'o-r')
title("Average")

figure(2)
plot(log10(N),theory_variance)
hold on
plot(log10(N),sample_variance,'o-r')
title("Variance")

figure(3)
plot(log10(N),ones([1 length(N)])*0.6)
hold on
plot(log10(N),prob,'o-r')
title("p < 60")

%% 1.3
clear all
close all
clc

N = [10E2 10E3 10E4 10E6 10E7];

h_1 = nan([1 length(N)]);
h_2 = nan([1 length(N)]);
i=1;
for N_i=N
x = rand(N_i,1);
y = rand(N_i,1);

z = x.*y;
h_1(i) = mean(z);
h_2(i) = mean(x)*mean(y);
i = i+1;
end

figure(1)
plot(log10(N),h_1,'.-')
hold on
plot(log10(N),h_2,'o-')
title("<z> = <x><y>")
legend(["mean(z)", "mean(x)*mean(y)"])