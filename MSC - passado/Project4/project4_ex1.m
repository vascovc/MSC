close all
clear variables
clc

H = 0.1;    %Magnetic Field
J = 1;      
N = 1000;   
m = 100000; 

T = 0.1:0.05:10;

M = zeros(m,1);
MT = zeros(length(T),1);
M100 = zeros(length(T),m/100-200);

S = ones(N,1);
l=1;


for k = 1:length(T)
    

    
    for i = 1:m

        p_center = randi(N);

        if p_center > 1 && p_center < N
            p_left = p_center - 1;
            p_right = p_center + 1;
        elseif p_center == 1
            p_left = N;
            p_right = p_center + 1;
        else
            p_left = p_center - 1;
            p_right = 1;
        end

        Spin        = S(p_center);
        SpinLeft    = S(p_left);
        SpinRight   = S(p_right);

        %% S Change
        E      = -J*Spin*(SpinLeft + SpinRight) - H*Spin;
        E_new   = -J*(-Spin)*(SpinLeft + SpinRight) - H*Spin;
        dE = E_new - E;

        if dE <= 0
            S(p_center) = -S(p_center);
        else  
            beta = 1/T(k);
            w = exp(-beta*dE);
            r = rand(1);
            if r <= w
                S(p_center) = -S(p_center);
            end
        end

        M(i) = (1/N)*sum(S);
    end
    
    % Average Magnetization

    for i = 100:100:m-100
        M100(k,l) = (1/i)*sum(M(1:i+99));
        l = l + 1;
    end
    
    MT(k) = mean(M100(k,:));
end


M_theo = zeros(length(T),1);
for j = 1:length(T)
    beta = 1/T(j);
    M_theo(j) = (sinh(beta*H)) / ( sqrt( (sinh(beta*H))^2 + exp(-4*beta*J)) );
end

%% Plotting
close all

figure
% set(gcf,'units','normalized','position',[0 0 0.8 0.4])
k_T=[1, 99, length(T)];
plot(M100(k_T(1),:),'b--','LineWidth', 1.2);
xlim([0 200]); xlabel('microstates'); ylabel('Magnetization');
grid on, hold on
plot(M100(k_T(2),:),'g--','LineWidth', 1.2);
xlim([0 200]);
k = length(T);
plot(M100(k_T(3),:),'r--','LineWidth', 1.2);
xlim([0 200]);
title('Magnetization vs microstates')
legend(['T = ' num2str(T(k_T(1)))], ['T = ' num2str(T(k_T(2)))], ['T = ' num2str(T(k_T(3)))])


figure
plot(T,M_theo,'k-');
hold on, grid on
plot(T,MT,'r.-');
xlabel('Temperature'); ylabel('Magnetization');
title('Magnetization vs Temperature')
legend('M_{theo}','M_{exp}')

figure
plot(log(T),log(M_theo),'k-');
hold on, grid on
plot(log(T),log(MT),'r.-');
xlabel('log_{10}Temperature'); ylabel('log_{10}Magnetization');
title('Magnetization vs Temperature (log-scale)')
legend('M_{theo}','M_{exp}')

s_MT=std(MT)
err=immse(MT,M_theo) %mean squared error

