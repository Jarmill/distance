classdef meas_wass_csp < handle
    %MEAS_WASS_CSP A Wasserstein joint probability measure for distance
    %estimation. Used to compute the distance between points on
    %trajectories and points on an unsafe set  
    %
    %The joint Wasserstein measure over (x, y) is split into n+1 measures
    %each in n variables (xk...xn, y1...yk) according to a correlative 
    %sparsity pattern. This assumes that the objective is separable into 
    %c(x, y) = sum_i ci(xi, yi). This is a very basic implementation of
    %CSP that black-boxes the structure of Y, and does not decompose
    %further based on constraint structure (e.g. Y is a box).    
    properties
        vars;       %variables in measures (@mpol)
        meas;       %cell array of measures (@meas)
        supp;
    end
    
    methods
        function obj = meas_wass_csp(vars, suffix, supp_X, supp_Xu)
            %UNTITLED2 Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
            
            obj.vars = struct('x', [], 'y', []);
            obj.vars.x = vars.x;
            obj.vars.y = vars.y;
            
            [obj.meas, obj.supp] = obj.meas_wass_def_csp(suffix, supp_X, supp_Xu);
        end
        
        function [meas_out, supp_out] = meas_wass_def_csp(obj,suffix, supp_X, supp_Xu)
            %MEAS_WASS_DEF_CSP define measures in this container
            %   Detailed explanation goes here
%             outputArg = obj.Property1 + inputArg;
            n = length(obj.vars.x);
            meas_out = cell(n, 1);
            supp_out = [];
            for i = 1:n
                name_curr = ['csp_', num2str(i), suffix];
 
                mpol(name_curr, n+1, 1);
                allvar_curr = eval(name_curr);

                
                %variables that are involved in the current measure   
                varcurr = struct('x',  allvar_curr(1:(end-i)), 'y', allvar_curr((end-i+1):end));             
                
                if i == 1
                    supp_curr = subs_vars(supp_X, obj.vars.x, varcurr.x);
                    supp_out = [supp_out; supp_curr];
                elseif i == n
                    supp_curr = subs_vars(supp_Xu, obj.vars.y, varcurr.y);
                    supp_out = [supp_out; supp_curr];
                else
                    supp_curr = 0;
                end                                
                
%                 meas_out{i} = obj.meas_def(suffix, varcurr, supp_curr);                
                meas_out{i} = meas_wass(varcurr, supp_curr);
            end

        end       
        
        %% get monomials 
        function mmmon_out = mom_monom_x(obj, dmin, dmax)
            %MOM_MONOM_X x-monomials (state) of wass csp measure
            %from degree dmin to dmax
            if nargin == 2
                dmax = dmin;
                dmin = 0;
            end
                        
            mmon_out = mmon(obj.meas{1}.vars.x, dmin, dmax);
            
            if isempty(obj.meas{1}.supp)
                %empty support: zero moments
                %may disable this
                mmon_out = zeros(size(mmon_out));
            end
            
            mmmon_out = mom(mmon_out);
        end    
        
        function mmmon_out = mom_monom_y(obj, dmin, dmax)
            %MMON monomials of variables of measure
            %from degree dmin to dmax
            if nargin == 2
                dmax = dmin;
                dmin = 0;
            end
                        
            mmon_out = mmon(obj.meas{end}.vars.y, dmin, dmax);
            
            if isempty(obj.meas{end}.supp)
                %empty support: zero moments
                %may disable this
                mmon_out = zeros(size(mmon_out));
            end
            
            mmmon_out = mom(mmon_out);
        end 
        
        %% constraint generator
        function con = overlap_con(obj, d)
           %form overlap constraints between cliques of the csp graph
           %this involves moment substitution, so dual variables will not
           %be recovered in the equality constraints. That is fine, because
           %this formulation will be equivalent to the CSP psatz
           %representation of sparse polynomials
           n = length(obj.vars.x);
           con = [];
           for i = 1:(n-1)
               
               %variables on overlap
               var_prev = [obj.meas{i}.vars.x(2:end); obj.meas{i}.vars.y]; 
               var_next = [obj.meas{i+1}.vars.x; obj.meas{i+1}.vars.y(1:end-1)]; 
               
               mom_prev = mom(mmon(var_prev, d));
               mom_next = mom(mmon(var_next, d));
               
%                con_curr = (mom_prev == mom_next);
               con_curr = (mom_prev - mom_next==0);
               
               con = [con; con_curr];               
           end
        end
        
        %% process the objective
        
        function mom_obj = mom_objective(obj, p, joint_vars)
            %MOM_OBJECTIVE find the csp moment representation of a
            %separable scalar term p with respect to the variables and
            %measures
            
            %joint_vars should be nothing
            
            p_public = get_representation(p, joint_vars);
            
            
            nterms = length(p_public.coef);
            n = length(obj.vars.x);
            %logical array determining which measure to use when writing
            %objective
            all_pow = reshape(p_public.pow, [nterms, n, 2]);
            screen_pow = any(all_pow, 3);
            
            %add up the terms in the objectives together
            mom_obj = 0;
            for i = 1:nterms                
                curr_screen = screen_pow(i, :);
                meas_active = obj.meas{curr_screen};
                curr_vars = [meas_active.vars.x(1); meas_active.vars.y(end)];
                curr_pow_x = all_pow(i, curr_screen, 1);
                curr_pow_y = all_pow(i, curr_screen, 2);
                
                curr_monom = prod(curr_vars.^[curr_pow_x; curr_pow_y]);
                
                mom_obj = mom_obj + mom(curr_monom*p_public.coef(i));
                
            end
            
        end
        
        %% post-processing
        function [optimal, mom_out, corner] = recover(obj, tol)
            %RECOVER if top corner of the wasserstein moment matrix is 
            % rank-1, then return approximate pair of points of closest 
            %  approach: (on trajectory, on unsafe set)
            
            if nargin < 2
                tol = 5e-4;
            end
            
            corner = cellfun(@(m) m.mmat_corner(), obj.meas, 'UniformOutput', false);
            rank_corner = cellfun(@(m) rank(m, tol), corner);
            n = length(obj.vars.x);
            if all(rank_corner == 1)
                %check the overlaps as well
                optimal = 1;
                mom_out.x = corner{1}(1+(1:n), 1);
                mom_out.y = corner{n}(2+(1:n), 1);
            else
                optimal = 0;                
                mom_out.x = [];                
                mom_out.y = [];
                
            end
%             corner = obj.mmat_corner();
%             
%             mass_curr= corner(1, 1);
%             if mass_curr < tol*1e-3
%                 %measure is empty
%                 %nobody home
%                 optimal = 1;                
%                 mom_out.x = [];                
%                 mom_out.y = [];
%                 corner = zeros(size(corner));
%             else
%                 rankM = rank(corner, tol);            
%                 optimal = (rankM == 1);
% 
%                 Nx = length(obj.vars.x);
%                 mom_out.x = corner(1+(1:Nx), 1);                
%                 mom_out.y = corner(1+Nx+(1:Nx), 1);                
% 
%             end
        end
        
                                
    end
    
end

