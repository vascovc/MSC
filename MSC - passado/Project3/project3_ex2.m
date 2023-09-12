clear all, clc

N=5e4; %nr of attempts
xc=-30;dt=10;
tf=50000;
D=1/4; %diffusion coefficient

t=0:tf;
F_c=zeros(1,length(t)); %first passage count
S_c=F_c; %survival count

for n=1:N
    r=zeros(length(t),2); %initialization the particle's positions 
    x=0;
    y=0;
    survival=1;
    for i=2:length(t)    % i is the number of steps from 1 to 10
        A=randi([1,4]);   % generates a random value drawn from the standard uniform distribution on the interval(0,1)
        mov_x=0;
        mov_y=0;
        if A==1
            mov_y=1;
        elseif A==2
            mov_x=1;
        elseif A==3
            mov_y=-1;
        elseif A==4
            mov_x=-1;
        end
        x=mov_x+x; 
        y=mov_y+y;
        if x==xc
            F_c(i)=F_c(i)+1;
            break
        else
            S_c(i)=S_c(i)+1;
        end
    end    
end

t_dt=1:dt:tf-2;
idx=1;
for i= 0:dt:length(t)-dt
    soma=0;
    soma2=0;
    for j=1:dt
        soma = soma + F_c(i+j);
        soma2 = soma2 + S_c(i+j);
    end
    N_f(idx)=soma;
    N_s(idx)=soma2;
    idx=idx+1;
end

%plotting results:

%first passage

F=N_f./(N*dt);
figure
subplot(121)
plot(t_dt,F,'.')
hold on, grid on


F_theo=abs(xc)./(sqrt(4*pi*D*t_dt.^3)) .* exp(-((xc)^2)./(4*D.*t_dt));

plot(t_dt,F_theo,'k-','LineWidth',1.5)
title('First passage Probability, F(t)')
legend('F','F_{theoretical}')
xlabel('t'),ylabel('F')
subplot(122)
plot(log10(t_dt),log10(F),'.')
hold on, grid on
plot(log10(t_dt),log10(F_theo),'k-','LineWidth',1.5)
%xlim([2 4])
title('First passage Probability, F(t) (log-log)')
legend('log_{10}(F)','log_{10}(F_{theoretical})')
xlabel('log_{10}t'),ylabel('log_{10}F')

S=(N_s./(N*dt));

figure

subplot(121)
plot(t_dt(2:end),S(2:end),'r--','LineWidth',3)
hold on, grid on
S_theo= erf(30./(2*sqrt(D.*t_dt)));
plot(t_dt,S_theo,'k-','LineWidth',1.5)
title('Survival Probability, S(t)')
legend('S','S_{theoretical}')

subplot(122)
plot(t_dt(2:end),S(2:end),'r--','LineWidth',3)
hold on, grid on
S_theo= erf(30./(2*sqrt(D.*t_dt)));
plot(t_dt,S_theo,'k-','LineWidth',1.5)
title('Survival Probability, S(t) [zoomed]')
legend('S','S_{theoretical}')
xlim([2500 2700])

