classdef meas_wass < meas_interface
    %MEAS_WASS A Wasserstein joint probability measure for distance
    %estimation. Used to compute the distance between points on
    %trajectories and points on an unsafe set    
    
    
    properties
%         Property1
        %vars: ('x', [], 'y', []):
        %   x: points on trajectories
        %   y: points on unsafe set
    end
    
    methods
        function obj = meas_wass(vars,supp)
            %MEAS_WASS Construct an instance of this class
            obj@meas_interface(vars, supp)
        end
        
        function mmmon_out = mom_monom_x(obj, dmin, dmax)
            %MMON monomials of variables of measure
            %from degree dmin to dmax
            if nargin == 2
                dmax = dmin;
                dmin = 0;
            end
                        
            mmon_out = mmon(obj.vars.x, dmin, dmax);
            
            if isempty(obj.supp)
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
                        
            mmon_out = mmon(obj.vars.y, dmin, dmax);
            
            if isempty(obj.supp)
                %empty support: zero moments
                %may disable this
                mmon_out = zeros(size(mmon_out));
            end
            
            mmmon_out = mom(mmon_out);
        end 
    
        function mom_obj = mom_objective(obj, p, joint_vars)
            %evaluate moments for location_distance/objective_con
            
            mom_obj = mom(obj.var_sub(joint_vars, p));
        end
        
        function mom_out= mom_sub(obj, vars_old, f_old)
            %MOM_PUSH pushforward moments v(f(x))
            
            mom_out = 0;
            for i = 1:length(obj.meas)
                curr_push = obj.var_sub(vars_old, f_old);
                mom_out = mom_out + mom(curr_push);
            end
        end  
        
        function [optimal, mom_out, corner] = recover(obj, tol)
            %RECOVER if top corner of the wasserstein moment matrix is 
            % rank-1, then return approximate pair of points of closest 
            %  approach: (on trajectory, on unsafe set)
            
            if nargin < 2
                tol = 5e-4;
            end
            
            corner = obj.mmat_corner();
            
            mass_curr= corner(1, 1);
            if mass_curr < tol*1e-3
                %measure is empty
                %nobody home
                optimal = 1;                
                mom_out.x = [];                
                mom_out.y = [];
                corner = zeros(size(corner));
            else
                rankM = rank(corner, tol);            
                optimal = (rankM == 1);

                Nx = length(obj.vars.x);
                mom_out.x = corner(1+(1:Nx), 1);                
                mom_out.y = corner(1+Nx+(1:Nx), 1);                

            end
        end
    end
end

