function R=restrict(n,N)
for i=1:n
    ind(2*i-1)=1+(i-1)*N;
    ind(2*i)=i*N;
end
R=eye(n*N);
R(ind,:)=[];