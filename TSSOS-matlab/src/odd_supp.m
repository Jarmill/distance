function osupp=odd_supp(n,supp)
%find the monomials in supp where any exponent in the monomial is odd
%candidates: [1,1],  [2, 1],   [4, 4, 3]
% i=1;
% lsize=size(supp);
% lo=lsize(1);
% indexb=1:lo;
% while lo>=i
%       bi=supp(indexb(i),:);
%       if length(bi(~mod(bi,2)))==n
%          indexb(i)=[];
%          lo=lo-1;
%       else
%           i=i+1;
%       end
% end
% osupp=supp(indexb,:);

%use vectorized matlab expressions to make this simpler

%find the set of supports where any exponent is odd
osupp = supp(any(mod(supp, 2), 2), :);

end