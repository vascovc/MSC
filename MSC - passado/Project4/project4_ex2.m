close all
clear variables
clc

%% Vars
N = 100;        %Number of spins
m = 5000;      %Number of microstates
H = 0%0.001;      %Magnetic Field
J = 1;          %Energy unit
T = 0.1:0.05:10;


M = zeros(length(T),m);
MT = zeros(length(T),1);
Beta = 1./T;


%% Cycle
for k = 1:length(T)
    
    System = ones(N,1);
    q=1;
    for i = 1:m
	
        pCenter = randi(N);

        Spin        = System(pCenter);
        SpinRestSum = sum(System) - Spin;
        
        % System Change
        E      = -(J/N)*Spin*SpinRestSum - H*Spin;
        ENew   = -(J/N)*(-Spin)*SpinRestSum  - H*Spin;
        EDelta(i) = ENew - E;

        if EDelta(i) <= 0
            System(pCenter) = -System(pCenter);
        else  
            w = exp(-Beta(k)*EDelta(i));
            r = rand(1);
            if r <= w
                System(pCenter) = -System(pCenter);
            end
        end

        M(k,i) = (1/N)*sum(System);
%         chi(k,i) = ( N/T(k) )*( mean(M(k,:).*M(k,:)) - mean(M(k,:))^2);
        

        if i>100
            if mod(i,100) == 0  % de 100 em 100 calcular <M>
                M100(q) = (1/i) * sum(M(k,:));  
%                 chi100(q) = (1/i) * sum(chi(k,:)); 
                q=q+1;
            end    
        end
        
    end
    clc

    AvgM(k) = (1/m)*sum(M100); %M(k,:) 
    chi(k) =  N/T(k)* mean( M(:).^2 - AvgM(k)^2)  ;
%     Avgchi(k) = (1/m)*sum(chi100);
end
plot(EDelta)
Beta=1./T;
M_theor=tanh(Beta*J.*AvgM+Beta*H);
chi_theor=Beta./(cosh(Beta*J.*AvgM+Beta*H).^2-Beta*J);

%% Plots
close all

%Dependency on temperature
figure
set(gcf,'units','normalized','position',[0 0 0.8 0.4])
subplot(121)
title('<M> vs Temperature')
hold on
plot(T,AvgM,'r.-');
plot(T,M_theor,'b-');
plot([J J],[-0.5 10],'k--','LineWidth',1.3)
ylim([-0.5e-3 20e-3])
text(1.2,17e-3,'T=T_c')
grid on
xlabel('Temperature'); ylabel('Magnetization');
legend('M_{exp}','M_{theo}')
subplot(122)
title('<M> vs Temperature: near T_c')
hold on
plot(T,AvgM,'r.-');
plot(T,M_theor,'b-','LineWidth',1.3);
plot([J J],[-0.5 10],'k--','LineWidth',1.3)
text(1.05,6.1e-3,'T=T_c')
ylim([0.5e-3 8e-3]),xlim([0 2])
legend('M_{exp}','M_{theo}')
grid on
xlabel('Temperature'); ylabel('Magnetization');
%%

% X=1/ |T-Tc| at T>Tc and X=1/ (2|Tc-T|) at T<Tc
% at T>Tc , and  at T a little bit below Tc

% figure
% set(gcf,'units','normalized','position',[0 0 0.8 0.4])
% subplot(121)
% plot(T,chi,'.-');
% hold on
% grid on
% title('Mag. Susceptibility vs Temperature')
% xlabel('Temperature'), ylabel('\chi')
% plot([J J],[0 1000000],'k--','LineWidth',1.3)
% text(1.2,740,'T=T_c')

figure
plot(T,chi,'.-');
ylim([0 20])
hold on
plot(T,chi_theor)
ylim([0 20])


hold on
grid on
title('Mag. Susceptibility vs Temperature: near T_c')
xlabel('Temperature'), ylabel('\chi')
text(1.05,65,'T=T_c')
plot(T,chi_theor,'-')
legend('\chi_{exp}','\chi_{theo}')
xlim([0.9 10])
plot([J J],[0 1000],'k--','LineWidth',1.3)
%%
%For individual temperature
figure
l=1;
id_k=[find(T==0.1) find(T==J) find(T==5) find(T==10)] %temperature of index 50
for k=id_k
    figure
for i = 100:100:m-100
    MShow(l) = mean(M(k,i:i+99)); 
    l = l + 1;
end
plot(MShow,'.-');
title(['Magnetization vs microstates (T= ' num2str(T(k)) ')'])
xlim([0 100]);
xlabel('microstates'); ylabel('Magnetization');
end


