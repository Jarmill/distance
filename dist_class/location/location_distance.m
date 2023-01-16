classdef location_distance < location & location_distance_interface
    %LOCATION_DISTANCE A location (space) of a dynamical system
    % Specifically designed to compare distances from points along
    % trajectories to points in an unsafe set
   
    
    methods
        function obj = location_distance(unsafe_supp, f, loc_id)
            %LOCATION_DISTANCE Construct an instance of this class
            %   Detailed explanation goes here
            if nargin < 3
                loc_id = [];   
            end
            
            %superclass constructor
            obj@location(unsafe_supp, f, 0, loc_id);
            obj@location_distance_interface(unsafe_supp, loc_id);
                        
        end     
        
        function [objective, cons_eq, cons_ineq, len_dual] = all_cons(obj, d)
            %ALL_CONS all constraints involving solely location
            %does not include sum of mass of initial measures
            %Output:
            %   cons_eq: equality constraints
            %   cons_ineq: inequality constraints (objective)
            
            %gather all constraints together
            liou = obj.liou_con(d);
            len_liou = length(liou);
            [abscont, len_abscont] = obj.abscont_box_con(d);
            
            [objective, cons_ineq, cons_eq_obj]= obj.objective_con();
            
            marg = obj.marg_wass_con(d);
            len_marg = length(marg);
            
            %package up the output
            len_dual = struct;
            len_dual.v = len_liou;
            len_dual.zeta = len_abscont;
            len_dual.w = len_marg;
            len_dual.beta = length(cons_ineq);
            
            %ensure this iss the correct sign
            cons_eq = [-liou; abscont; marg]==0; 
            cons_eq = [cons_eq; cons_eq_obj];
        end
        
        function [objective, cons_ineq, cons_eq_obj] = objective_con(obj)
            [objective, cons_ineq, cons_eq_obj] = ...
                objective_con@location_distance_interface(obj);
        end
        
        function [len_out] = len_eq_cons(obj)
            %LEN_EQ_CONS Number of equality constraints strictly in this
            %location 
            len_out = obj.len_dual.v + sum(obj.len_dual.zeta) + obj.len_dual.w;
        end
        
        function obj = dual_process(obj, d, rec_eq, rec_ineq, gamma)
             %DUAL_PROCESS turn the dual variables from solution into 
             %polynomials and interpretable quantities
             %
             %Input:
             %  d:          2*order, order of auxiliary polynomials
             %  rec_eq:     dual variables from equality constraints
             %  rec_ineq:   dual variables from inequality constraints
             %  gamma:      objective value (as a dual variable)
             %TODO: fix this so that boxes can be used
             
             %numeric quantities
             obj.dual.solved = 1;
             
%              obj.dual.beta = rec_ineq;
             obj.dual.gamma = gamma;
             
             beta_count = 0;      
             if obj.len_dual.beta == 0
                 obj.dual.beta = 1;
             else
                 Nterms = length(obj.dist);
                 obj.dual.beta = cell(Nterms, 1);
                 for i =1:Nterms
                     len_dist = length(obj.dist{i});
                     obj.dual.beta{i} = rec_ineq(beta_count + (1:len_dist));
                     beta_count = beta_count + len_dist;
                 end
             end
             
             %process the polynomials
             
             %auxiliary function v
             v_coeff = rec_eq(1:obj.len_dual.v);
             monom = mmon(obj.get_vars_end(), 0, d);
             obj.dual.v = v_coeff'*monom;
             
             count_zeta = obj.len_dual.v;
             
             %iterate through all subsystems
             Nb = length(obj.vars.b);
             monom_all = mmon(obj.get_vars(), 0, d);
             
             %TODO: confirm that all abscont relations have the same length
             len_monom_all = length(monom_all);
             for i = 1:length(obj.sys)       
                 
                 %untangle the box zeta functions
                 zeta = [];
                 for j = 1:Nb
                     zeta_coeff = rec_eq(count_zeta + (1:len_monom_all));
                     zeta_curr = zeta_coeff'*monom_all;
                     zeta = [zeta; zeta_curr];   
                     
                     count_zeta = count_zeta + len_monom_all;
                 end
                 

                 %ship off dual variables for processing in subsystem 
                 obj.sys{i} = obj.sys{i}.dual_process(obj.dual.v, zeta);       
                 
                 
                 %figure out nonnegativity requirement
             end
             
             %take care of the marginal constraint w(x)
             w_coeff = rec_eq(count_zeta + (1:obj.len_dual.w));
%              x_monom = 
             x_monom = mmon(obj.vars.x, 0, d);
             obj.dual.w = w_coeff'*x_monom;
            
             %nonnegativity of location (not subsystems)
            %initial measure
            if isempty(obj.init)
                nn_init = 0;
            else
                nn_init = -obj.dual.gamma + obj.dual.v;
            end
            
            %terminal measure
            if isempty(obj.term)
                nn_term = 0;
            else
%                 nn_term = obj.dual.v - obj.dual.beta'*obj.objective;
                nn_term = obj.dual.w - obj.dual.v;
            end
            
            if obj.len_dual.beta == 0
               nn_wass = obj.dist - obj.dual.w; 
            else
                nn_wass = -obj.dual.w;
                for i = 1:length(obj.dist)
                    nn_wass = nn_wass + obj.dual.beta{i}'*obj.dist{i};
                end
            end
            
            %process the wasserstein dual
            obj.dual.nn = [nn_init; nn_term; nn_wass];
             
        end      
        
        %% support 
        function supp_con_out = supp_con(obj)
            supp_con_out = supp_con@location(obj);
            wass_con = supp_con@location_distance_interface(obj);
            
            supp_con_out = [supp_con_out; wass_con];
        end
        
        function [optimal, mom_out, corner] = recover(obj, tol)
            %RECOVER if top corner of the moment matrix is rank-1, then
            %return approximate optimizer
            if nargin < 2
                tol = 5e-4;
            end

%             
            [optimal, mom_out, corner] = recover@location_interface(obj, tol);
            [optim_wass, mom_wass, corner_wass] = recover@location_distance_interface(obj, tol);
            
            optimal = optimal && optim_wass;
            corner.wass = corner_wass;
            mom_out.y = mom_wass.y;

        end
            
        
        
    end
end

