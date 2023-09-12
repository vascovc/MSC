%%% Nuno Monteiro (79907)

%%% Modelaçao de Sistemas Complexos
%% task 1

clear all, close all, clc

%%%%% task 1.2

M=100
idx=1;
%analytical results:
mean_x_an=(M+1)/2;
var_x_an=(M^2-1)/12;

for N = [1e3,1e4,1e6]
  x=randi(M,[1 N]); %generating N random numbers with unif distr.
  mean_x(idx)=sum(x)/N;
  var_x(idx)=sum((x-mean_x(idx)).^2)/N;
  
  %calculation of relative errors (%)
  re_mean_x(idx)=abs(mean_x(idx)-mean_x_an)/mean_x_an*100;
  re_var_x(idx)=abs(var_x(idx)-var_x_an)/var_x_an*100;
  
  idx=idx+1;
end
mean_x
var_x
re_mean_x
re_var_x

%%% task 1.3

clear all, clc

idx=1;
m_z_1=zeros(3,1);
m_z_2=m_z_1;
re=m_z_1;
for N = [1e2,1e4,1e6]
  x=zeros(N,1);
  y=x;
  z=x;
    for i=1:N
      x(i)=rand(1);
      y(i)=rand(1);
      z(i)=x(i)*y(i);
    end
  m_x=mean(x);
  m_y=mean(y);
  m_z_1(idx)=mean(z);
  m_z_2(idx)=m_x*m_y;
  re(idx)=(abs(m_z_1(idx)-m_z_2(idx))/m_z_1(idx))*100;
  idx=idx+1;
end
m_z_1
m_z_2
re
% 
%% task 2

clear all, close all, clc

%%%%%%%%%%%%%%%%%
%%% task 2.1
N=1e5
x=rand(N,1);
dx = 0.005;

figure('units','normalized','outerposition',[0 0 0.7 0.6])
subplot(121)
hist(x,30),title('(a) - Histogram of 10⁵ x=rand(1) values')
ylabel('Frequency'),xlabel('x'),ylim([0 4000])

[M,x_m,bins] = histcounts(x,1/dx);
P=M./(N*dx);
normalize=sum(P*dx) %checking normalization (=1)

subplot(122)
plot(P,'k.-'),ylim([0 1.2])
ylabel("P(x)"),xlabel("x"),title("(b) - PDF of x (dx=1/30)")
print('./latex/Figures/fig1_task2','-dpng')

%%%%%%%%%%%%%%%%%%%
%%% task 2.2

m_x=mean(x);
m_x
v_x=var(x);
v_x

%accuracies: rel. error in %
re_m_x=100*abs(m_x-0.5)/0.5;
re_m_x 
re_v_x=100*abs(v_x-1/12)/(1/12);
re_v_x

%%%%%%%%%%%%%%%%%%
%%% task 2.3

N = [100,1e3,1e4,1e5,1e6];


%anallytical results:

for n=1:length(N)
    x = sqrt(rand(N(n),1));  

    [M,X_m,bins] = histcounts(x,(1/dx));

    P_n(:,n) = M./(N(n)*dx);
    for k = 0:1/dx-1
        x_k(k+1)=(k+0.5)*dx;  %x axis
    end
    
    %theoretical results
    g=2*sqrt(x_k);
    t_mean_x=mean(g)
    t_var_x=var(g);

    normalize(n)=sum(P_n(:,n).*dx);
    
    %numerical results:
    mean_x(n)=2*sum((x_k.*P_n(:,n)').*dx); %why adjust to multiplication to 2?
    std_x(n)=std(x)
    if n==1
        figure('units','normalized','outerposition',[0 0 0.4 1])
    end
    while 1
        if n==1
            subplot(3,1,1)
        elseif n==2
           subplot(3,1,2) 
        elseif n==5
           subplot(3,1,3) 
        else
            break
        end
        plot(x_k, P_n(:,n), 'k.-'); hold on
        plot(sqrt(x_k), 2*sqrt(x_k),'k-')
        title(['Probability density Function, N=',num2str(N(n))])
        legend('Numerical PDF','Anallytical PDF');
        xlabel('x'); ylabel('P');
        break
    end
end
print('./latex/Figures/fig2_task2','-dpng')

figure
plot([-1 7],[t_mean_x t_mean_x],'LineWidth',3),hold on
plot(log10(N),mean_x,'ko')
xlim([1.5 6.5]),ylim([1.1 1.4])
title('Averaged values of x')
legend('Theoretical mean','Numerical mean')
xlabel('log_{10}(N)'),ylabel('<x>')
print('./latex/Figures/fig3_task2','-dpng')

acc_mean=abs(t_mean_x-mean_x)/t_mean_x * 100

%% task 3

clear all, close all, clc
dy = 0.005;
m=1e6
cols=['r','k','b'];
n=[10,100,1000];
for j=1:3
    nn=n(j);
    dy=0.005;
    Y=zeros(m,1);
    for i=1:m
        z=rand(1,nn);
        Y(i)=sum(z)/nn;
    end
    % Find the bin placement based on the Y values
    [M,Y_m,bins] = histcounts(Y,(1/dy)-1);

    P_yk=zeros(1/dy-1,1)
    for k = 1:1/dy-1
        P_yk(k)=M(k)/m;
        y_k(k)=(k+0.5)*dy;
    end
    normalize=sum(P_yk) %checking normalization condition
    
    if j==1
        figure('units','normalized','outerposition',[0 0 0.7 0.6])
    end
    subplot(121)
    plot(y_k,P_yk,[cols(j),'-'])
    hold on
    mean_Y(j)=mean(Y_m);
    var_Y(j)=var(Y_m);
    sd_Y(j)=(std(Y_m)^2)/nn;
    ylabel("P"),xlabel("Y"),title("Probability Density of x"),legend(['n=' num2str(n(1))],['n=' num2str(n(2))],['n=' num2str(n(3))])
    ylim([0 0.025])
    
end
subplot(122)
plot(log10([10,100,1000]),var_Y,'ko-')
hold on, grid on
plot(log10([10,100,1000]),sd_Y,'k.-','MarkerSize',15)
legend('var(Y)','\sigma^2/n')
xlabel('log_{10}(n)')
ylabel('Y')
xlim([0.7 3]),ylim([-0.02 0.1])
title('Variance and \sigma^2/n of Y')
print(['./latex/Figures/fig1_task3'],'-dpng')
mean()
%% task 4

clear all, close all, clc

M=9;N=21; %boxes, balls
j=1;
P_n=[];
K=[1e3,1e4,1e6] %nr of trials
p=1/M;
k=1;
N_tr=zeros(N+1,length(K));
P_n=N_tr;
markers=['*','o','s']
figure('units','normalized','outerposition',[0 0 0.7 0.7])
for k=1:length(K)
    n=zeros(K(k),1);
    for j=1:K(k)
        x=randi(M,[1 N]);
        n(j)=sum(x==3);
        N_tr(n(j)+1,k)=N_tr(n(j)+1,k)+1;
    end
    for i=1:N+1
        P_n(i,k)=N_tr(i,k)/K(k);
        i=i-1;
        Bin(i+1,k)=comb(i,N)*(p^i)*(1-p)^(N-i)
        Poi(i+1,k)=exp(-N*p)*((N*p).^i)/factorial(i);
        G(i+1,k)=(1/sqrt(2*pi*N*p))*exp(-((i-N*p)^2)/(2*N*p));
        i=i+1;
    end
    subplot(221), hold on
    plot(0:N,P_n(:,k),['k-' markers(k)]),title(["Observed P(n) in box nr. 3"])
    xlabel('n'), ylabel('P_n')
    legend(['K=' num2str(K(1))],['K=' num2str(K(2))],['K=' num2str(K(3))])
    subplot(222), hold on
    plot(0:N,P_n(:,k),['k-' markers(k)]),title(["[Expanded] Observed P(n) in box nr. 3"])
    xlim([0 4]),ylim([0.15 0.35])
    xlabel('n'), ylabel('P_n')
    legend(['K=' num2str(K(1))],['K=' num2str(K(2))],['K=' num2str(K(3))])
    if k==3
        subplot(223), hold on
        plot(0:N,Poi(:,k),'k-*')
        plot(0:N,G(:,k),'k.-')
        plot(0:N,Bin(:,k),'ko-')
        plot(0:N,P_n(:,k),'r-'),title(["Observed P(n) in box nr. 3"])
        xlabel('n'), ylabel('P_n')
        legend("Poisson","Gaussian","Binomial","Obs. vals.")
        title(["Theoretical P(n) vs P(n), K=10^6"])
        subplot(224), hold on
        plot(0:N,Poi(:,k),'k-*')
        plot(0:N,G(:,k),'k.-')
        plot(0:N,Bin(:,k),'ko-') 
        plot(0:N,P_n(:,k),'r-')
        xlim([0 4]),ylim([0.15 0.35])
        xlabel('n'), ylabel('P_n')
        legend("Poisson","Gaussian","Binomial","Obs. vals.")
        title(["[Expanded] Theoretical P(n) vs P(n), K=10^6"])
    end

end
print(['./latex/Figures/fig1_task4'],'-dpng')
%correlations
mean(corr(P_n(:,3),G(:,3)))
mean(corr(P_n(:,3),Poi(:,3)))
mean(corr(P_n(:,3),Bin(:,3)))