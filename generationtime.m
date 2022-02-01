
% This function computes the three measures A, T and mu and the population
% growth rate. The inputs are b=(b_1,...,b_n), s=(s_1,...,s_n) and
% S=(S_1,...,S_n). Here b_j is the fecundity at age j, s_j is the survival
% probability from age j-1 to age j and S_j=s_1*s_2*...s_j is the survival
% probability from birth to age j. Note that you can also compute s from S:
% s_1=S_1 ans s_j=S_j/S_{j-1} for j=2,...,n. The number n is the upper
% limit for individual life length in the population under study. (The
% reason for requiring input of both s and S, even though one can be
% computed from the other, is to avoid the need for one function when you
% want to input s and  another when you want to input S.)


function [A,T,mu,lambda] = generationtime(b,S,s)

n=length(b);

f=b.*s; % computes fertility from fecundity and survival rate.

M=[diag(s(1:n-1)), zeros(n-1,1)];
L=[f;M]; %this row and the previous compute the Leslie matrix
lambda=eig(L);
lambda=lambda(1);
p(1)=f(1);
p(2:n)=S(1:n-1).*f(2:n); % computes the p-vector
R0=sum(p); 

j=1:n;
v1=j.*p;
v2=lambda.^(-j);

A=v1*v2';
mu=j*p'/R0;
T=log(R0)/log(lambda);

