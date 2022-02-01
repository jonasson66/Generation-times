
% This script allows you to pick fecundity function, survival rate and
% corresponding according to the standard function given in this paper. 
% The choice is made by uncommenting the corresponding rows.
% The script then outputs the generation time measures A, T and mu and 
% the population growth rate lambda. 
% These are here given the names B, U, M and L respectively in order for
% looping over a parameter of choice to produce vectors A, T, mu, lambda,
% where each coefficient corresponds the corresponding value of the
% parameter that is varied.

j=1:n;

b=gamma*ones(1,n);
%b(1:4)=zeros(1,4);

%b=gamma*exp(-kappa*j);

%b=gamma*(1-exp(-kappa*j)).^3;

%b=zeros(1,n);, b(1)=c;, b(20)=c;


%s=exp(-alpha)*exp(-beta*exp(-rho*j)*(exp(rho)-1));
%S=exp(-alpha*j).*exp(-beta*(1-exp(-rho*j)));
 
%s=exp(-alpha*exp(beta*j)*(1-exp(-beta)));
%S=exp(-alpha*(exp(beta*j)-1));

s=exp(-alpha)*ones(1,n);
S=exp(-alpha*j);

[B,U,M,L]=generationtime(b,S,s);
