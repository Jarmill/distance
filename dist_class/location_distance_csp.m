classdef location_distance_csp < location_distance
    %LOCATION_DISTANCE_CSP A location (space) of a dynamical system
    % Specifically designed to compare distances from points along
    % trajectories to points in an unsafe set
    %
    % Uses a correlative sparsity pattern to decompose the joint constraint
    % over (x, y) into n+1 measures over (x(i:n), y(1:i))
    
    properties
%         Property1
    end
    
    methods
        function obj = location_distance_csp(unsafe_supp, f, id)
            %LOCATION_DISTANCE_CSP Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3
                id = [];   
            end
            
            obj@location_distance(unsafe_supp, f, id);
        end
        
        
        %create the distance measures        

        function meas_new = meas_wass_def(obj, suffix, supp_X, supp_Xu)
            meas_new = meas_wass_csp(obj.vars, suffix, supp_X, supp_Xu);
        end
        
        %constraint generation
        function [cons, len_w] = marg_cons(obj, d)
            [cons, len_w] = marg_cons@location_distance(obj, d);
            
            cons = [cons; obj.overlap_cons(d)];
        end

        function cons = overlap_cons(obj, d)
            cons = [];
            for i = 1:length(obj.wass)
                cons = [cons; obj.wass{i}.overlap_con(d)];
            end
        end

        function [objective, cons_eq, cons_ineq, len_dual] = all_cons(obj, d)
            [objective, cons_eq, cons_ineq, len_dual] = all_cons@location_distance(obj, d);
            
            %include the csp clique overlap constraints
            cons_eq = [cons_eq; obj.overlap_cons(d)];
        end
    end
end

