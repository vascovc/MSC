%% exercise 1.1

clear all, close all, clc

tf=50;
t=0:tf;
p = 0.5; %prob to jump on the left or on the right (p=q)
pos(1)=0; % 1-d position starts on the origin
figure
for n=1:3
for i=1:length(t)-1    % i is the number of steps from 1 to 10
    R=rand;   % generates a random value drawn from the standard uniform distribution on the interval(0,1)
    if R<p
        S=-1;
    elseif R>p
        S=1;
    end
    pos(i+1)=S+pos(i); % new position of this random walk
end
plot(pos,t,'o-'),grid on
hold on
end
xlim([-20 20])
xlabel('Position')
ylabel('t (steps)')
title('1-D Random Walks (t_f=50)')
legend('trajectory 1','trajectory 2','trajectory 3')
print('-djpeg90','./latex/Graphics/randwalk.jpg')


%% exercise 1.2
clear all, close all, clc
tf_matrix=[40,41;400,401;4000,4001]
p = 0.5; %prob to jump on the left or on the right (p=q)
N=50000; % nr of trajectories
x=-500:500;
P_mean=zeros(length(x),3);
for o=1:3
P=zeros(length(x),2);
tf_list=tf_matrix(o,:); % pairs considered: 40,41;400,401;4000,4001
for n=1:N
    for j=1:length(tf_list)
        tf=tf_list(j);
        t=0:tf;
        pos=zeros(length(t),1);
        for i=1:length(t)-1    % i is the number of steps from 1 to 10
            R=rand;   % generates a random value drawn from the standard uniform distribution on the interval(0,1)
            if R<p
                S=-1;
            elseif R>p
                S=1;
            end
            pos(i+1)=S+pos(i); % new position of this random walk
            if i==length(t)-1
                xi=pos(i+1);
                idx=xi+abs(min(x))+1;
                P(idx,j)=P(idx,j)+1;
            end
        end
    end
end

for i=1:length(P)
    P_mean(i,o)=(P(i,1)+P(i,2))/2;

end
P_mean(:,o)=P_mean(:,o)./N;

end
figure('units','normalized','outerposition',[0 0 0.5 0.6])
plot(x,P_mean(:,1),'r',x,P_mean(:,2),'g',x,P_mean(:,3),'b','LineWidth',1.5)
title('<P(x,t)> of 1-D Random Walks')
legend(['t = ' num2str(tf_matrix(1,1)) ', ' num2str(tf_matrix(1,2))] ...
        ,['t = ' num2str(tf_matrix(2,1)) ', ' num2str(tf_matrix(2,2))] ...
        ,['t = ' num2str(tf_matrix(3,1)) ', ' num2str(tf_matrix(3,2))])
xlabel('x')
ylabel('Probabilities')
xlim([-200 200])
grid on
print('-djpeg90','./latex/Graphics/p_3times.jpg')

tf=(tf_list(1)+tf_list(2))/2;
P_theoretical=1/sqrt(2*pi*tf).*exp(-(x.^2)/(2.*tf));
figure('units','normalized','outerposition',[0 0 0.5 0.6])
grid on, hold on
xlim([-200,200])
plot(x,P_mean(:,3),x,P_theoretical,'-','LineWidth',1.5)
title('Theoretical P(x,t) and <P(x,t)> (t=4000,4001)')
legend('<P(x,t)>','P_{t}')
xlabel('x')
ylabel('Probabilities')
print('-djpeg90','./latex/Graphics/p_comparison.jpg')
%mean squared errors
err1 = immse(P_theoretical',P_mean(:,1)) 
err2 = immse(P_theoretical',P_mean(:,2)) 
err3 = immse(P_theoretical',P_mean(:,3)) 

%% exercise 2

clear all, close all, clc

delta=0.015
p=0.5-delta
q=0.5+delta
tf_matrix=[40,41;400,401;4000,4001]
N=50000; % nr of trajectories
x=-500:500;
P_mean=zeros(length(x),3);
for o=1:3
P=zeros(length(x),2);

tf_list=tf_matrix(o,:); % pairs considered: 40,41;400,401;4000,4001
for n=1:N
    for j=1:length(tf_list)
        tf=tf_list(j);
        t=0:tf;
        pos=zeros(length(t),1);
        for i=1:length(t)-1    % i is the number of steps from 1 to 10
            R=rand;   % generates a random value drawn from the standard uniform distribution on the interval(0,1)
            if R<p
                S=-1;
            elseif R>p
                S=1;
            end
            pos(i+1)=S+pos(i); % new position of this random walk
            if i==length(t)-1
                xi=pos(i+1);
                idx=xi+abs(min(x))+1;
                P(idx,j)=P(idx,j)+1;
            end
        end
    end
end

for i=1:length(P)
    P_mean(i,o)=(P(i,1)+P(i,2))/2;

end
P_mean(:,o)=P_mean(:,o)./N;

end
figure('units','normalized','outerposition',[0 0 0.5 0.6])
plot(x,P_mean(:,1),'r',x,P_mean(:,2),'g',x,P_mean(:,3),'b','LineWidth',1.5)
title('<P(x,t)> of 1-D Asymmetric Random Walks')
legend(['t = ' num2str(tf_matrix(1,1)) ', ' num2str(tf_matrix(1,2))] ...
        ,['t = ' num2str(tf_matrix(2,1)) ', ' num2str(tf_matrix(2,2))] ...
        ,['t = ' num2str(tf_matrix(3,1)) ', ' num2str(tf_matrix(3,2))])
xlabel('x')
ylabel('Probabilities')
xlim([-200 350])
grid on
print('-djpeg90','./latex/Graphics/p2_3times.jpg')

tf=(tf_list(1)+tf_list(2))/2;
P_theoretical=1/sqrt(2*pi*tf).*exp(-((x-2*tf*delta).^2)./(2.*tf));
figure('units','normalized','outerposition',[0 0 0.5 0.6])
grid on, hold on
xlim([-200,350])
plot(x,P_mean(:,3),x,P_theoretical,'-','LineWidth',1.5)
title('Theoretical P(x,t) and <P(x,t)> (t=4000,4001)')
legend('<P(x,t)>','P_{t}')
xlabel('x')
ylabel('Probabilities')
print('-djpeg90','./latex/Graphics/p2_comparison.jpg')

% sum P =1 (checking 1st property in step 2 of the algorithm)
s1=sum(P_mean(:,1))
s2=sum(P_mean(:,2))
s3=sum(P_mean(:,3))

%checking 2nd property in step 2 of the algorithm
s1_2=sum(P_mean(:,1)).*(tf_matrix(1,1)+tf_matrix(1,2))*delta
s2_2=sum(P_mean(:,2).*(tf_matrix(2,1)+tf_matrix(2,2))*delta)
s3_2=sum(P_mean(:,3).*(tf_matrix(3,1)+tf_matrix(3,2))*delta)

%checking 3rd property in step 2 of the algorithm
s1_2=sum((P_mean(:,1)-(tf_matrix(1,1)+tf_matrix(1,2))*delta).^2)
s2_2=sum((P_mean(:,2)-(tf_matrix(2,1)+tf_matrix(2,2))*delta).^2)
s3_2=sum((P_mean(:,3)-(tf_matrix(3,1)+tf_matrix(3,2))*delta).^2)

%mean squared errors
err1 = immse(P_theoretical',P_mean(:,1)) 
err2 = immse(P_theoretical',P_mean(:,2)) 
err3 = immse(P_theoretical',P_mean(:,3)) 

%maximum values and comparing to the theory
max1=max(P_mean(:,1))
max2=max(P_mean(:,2))
max3=max(P_mean(:,3))
max1t=1/sqrt(2*pi*40.5)
max2t=1/sqrt(2*pi*400.5)
max3t=1/sqrt(2*pi*4000.5)
