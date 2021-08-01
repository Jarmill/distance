%from https://github.com/wangjie212/TSSOS/tree/matlab
%using symengine
syms x1 x2;
f=x1^4+x2^4-x1*x2;
g_1=1-x1^2-2*x2^2;
g_2 = 1 - x1*x2^2;
n=2;
m=2;
d=3; % the order of Lasserre's hierarchy
dg=[2, 3]; % the degree vector of {g_j}
[opt,data,status]=blockpop_cons_first(f,[g_1, g_2],n,m,d,dg,'mosek')