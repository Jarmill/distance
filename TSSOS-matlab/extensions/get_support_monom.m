function [data] = get_support_monom(p, X, vars, d)
%GET_SUPPORT_MONOM Find information about the support structure of a tssos
%psatz program. This is written for Yalmip symvar, future work may involve
%the symengine in syms
%
%Input:
%   p:      A polynomial that should be nonnegative on a region
%   X:      Region of nonnegativity (struct with fields ineq and eq)
%   vars:   variables of polynomial
%   d:      degree of polynomial (after multiplication by constraint)
%
%Output:
%  data: information about the support structure of (p, X)
%        each attribute is a cell array in terms of [p, g1>=0, g2>=0...]
%        (polynomial on index 1 and inequality constraints on index 2+)
%           fbasis: monomials on row of moment matrix (orders up to d)
%           gbasis: monomials on row of localizing matrix (orders up to d - deg(gk))
%           lt:     cardinality of support
%           coe:    coefficients of monomials
%           ssupp:  monomial powers present
%           supp:   all monomial powers (union over each ssupp)

%some code borrowed from src/blockpop_cons_first.m by Jie Wang

%% initializing fields
n = length(vars);   %number of variables
m = length(X.ineq); %number of inequality constraints

%all polynomials used (ignoring equalities)
pop = [p; X.ineq];

coe=cell(1,m+1);    %coefficients
ssupp=cell(1,m+1);  %monomial support
lt=zeros(1,m+1);    %size of support
deg = zeros(1, m+1);

%% process the polynomial p >= 0

for k = 1:m+1
    [ssupp{k}, coe{k}] = getexponentbase(pop(k), vars);
    lt(k) = length(coe{k});   
    deg(k) = max(sum(ssupp{k}, 2));
end

ssupp{1}=unique(ssupp{1},'rows', 'stable');
lt(1) = size(ssupp{1}, 1);

%% generate bases

%newton polytopes are not used in constrained optimization
%the full basis is used instead

fbasis=get_basis(n,d);
flb=nchoosek(n+d,d);
glb=zeros(1,m);
gbasis=cell(1,m);
for k=1:m
    glb(k)=nchoosek(n+d-ceil(deg(k+1)/2),d-ceil(deg(k+1)/2));
    gbasis{k}=get_basis(n,d-ceil(deg(k+1)/2));
end


%% package up to output
data = struct;
data.fbasis=fbasis;
data.gbasis=gbasis;
data.lt=lt;
data.coe=coe;
data.ssupp=ssupp;

% data.fsizes=fsizes;
% data.fsupp=fsupp;
% data.gsupp=gsupp;





end

