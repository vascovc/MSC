clear
clc
N = 10000;

for c=[50,100,500]
disp(['c=',num2str(c)])
L = N/2*c;
A = create_adjacency_matrix(N,L);

tic
ntr = trace(A^3);
ntr = ntr*1/6;
toc
disp(num2str(ntr))

tic
ntr=0;
for i=1:N
    for j=i:N
        if A(i,j)
            for k=j:N
                if A(j,k)==1 && A(k,i)==1
                    ntr = ntr+1;
                end
            end
        end
    end
end
toc
disp(num2str(ntr))
end