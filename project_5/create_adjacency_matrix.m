function A = create_adjacency_matrix(N,L)
A = zeros(N);
counter = 0;
while counter<L
    i = round(1 + (N-1).*rand);
    j = round(1 + (N-1).*rand);
    if A(i,j)~=1 && i~=j
        A(i,j) = 1;
        A(j,i) = 1;
        counter = counter+1;
    end
end
end