function [cliques] = find_csp(p, X, vars)
%FIND_CSP find correlative sparsity pattern associated with p>=0 over
%X.ineq >= 0 and X.eq == 0 

%find monomials of vars in p

[exp_p, base_p] = getexponentbase(p, vars);

G_p = (exp_p'*exp_p);

%TODO: reduction according to X.eq == 0

%find variables that are together in constraints X.ineq >= 0
exp_ineq = zeros(length(X.ineq), length(vars));
for i = 1:length(X.ineq)
    [exp_curr, base_curr] = getexponentbase(X.ineq(i), vars);
    exp_ineq(i, :) = any(exp_curr, 1);
%     exp_ineq = exp_ineq
end
G_ineq = (exp_ineq'*exp_ineq);

%this is the correlative sparsity graph
G = G_p + G_ineq;
G(G>0) = 1;
%find cliques of chordal extension of G
cliques = cliquesFromSpMatD(G);

%TODO: find cliques through maxCardinalitySearch.m if data is G
%chordal. The cholesky chordal completion method can fail

%TODO: assign constraints to cliques




end

