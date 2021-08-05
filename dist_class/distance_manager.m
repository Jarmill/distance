classdef distance_manager < manager_interface
    %DISTANCE_MANAGER Summary of this class goes here
    %   Detailed explanation goes here    
    methods
        function obj = distance_manager(unsafe_supp, f)
            %DISTANCE_MANAGER Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;

            %set up the location
            if unsafe_supp.CSP
                loc_curr = location_distance_csp(unsafe_supp, f, []);
            else
                loc_curr = location_distance(unsafe_supp, f, []);
            end
            obj@manager_interface(loc_curr);
                        
        end
        
  %% Formulating and solving program
        function sol = solve(obj, objective, mom_con,supp_con)
            %SOLVE formulate and solve peak estimation program from
            %constraints and objective    
            
            %TODO: incorporate minquery into maximin (minimax) formulation

            mset('yalmip',true);
            %make sure the solution is precise
            mset(obj.sdp_settings);
            % https://docs.mosek.com/9.2/pythonfusion/parameters.html
            
            

            P = msdp(min(objective), mom_con, supp_con);

            sol = struct;
            [sol.status,sol.obj_rec, ~,sol.dual_rec]= msol(P);        
        end  
        
        function [objective, mom_con, supp_con, len_dual] = cons(obj,d, Tmax)
            %PEAK_CONS formulate support and measure constraints for peak
            %program at degree d
            %Input:
            %   d:      Monomials involved in relaxation (2*order)
            %   Tmax:   Maximum time (only when time-independent)
            %
            %Output:
            %   objective:  target to maximize  (@mom)
            %   mom_con:    moment constraints  (@momcon)
            %   supp_con:   support constraints (@supcon)
          

            supp_con = obj.loc.supp_con();       %support constraint                 
                        
            %gather all constraints in each location
            %(loop over locations)
            [objective, cons_eq, cons_ineq, len_dual] = obj.loc.all_cons(d);
            %finalize moment constraints
            
            %mass of initial measure sums to one
            mass_init_con = (obj.loc.mass_init() - 1 == 0);
            
            if obj.loc.supp.TIME_INDEP
                mass_occ_con = (obj.loc.mass_occ() <= Tmax);
            else
                mass_occ_con = [];
            end

            %time independent: mass of sum of occupation measures are less
            %than Tmax. Implement/fix this?
            %TODO: get this done
                
%             len_liou = length(liou_con);

            mom_con = [cons_eq; cons_ineq; mass_init_con; mass_occ_con];
            

        end  
        
        function obj = dual_process(obj, d, dual_rec, len_dual)
            %DUAL_PROCESS dispatch the dual variables from solution to
            %locations and measures, turn the variables into nonnegative
            %functions along trajectories
            
            %v: coefficients of liouville
            %beta: coefficients of cost (if able)
            %alpha: dual of zeno gaps
            
            rec_eq = dual_rec{1};
            rec_ineq = dual_rec{2};
            
            %liouville
            %time independent
%             v_coeff = rec_eq(1:len_liou);
                        
            gamma = rec_eq(end);
            
            %counters for constraints (for multiple locations)
            count_eq = 0;
            count_ineq = 0;
            
            %index out current dual variable coefficients
            obj.loc.len_dual = len_dual;
            len_eq_curr = obj.loc.len_eq_cons();
            len_ineq_curr = sum(obj.loc.len_dual.beta);
            
            rec_eq_curr = rec_eq(count_eq + (1:len_eq_curr));
            rec_ineq_curr = rec_ineq(count_ineq + (1:len_ineq_curr));
            
            obj.loc = obj.loc.dual_process(d, rec_eq_curr, rec_ineq_curr, gamma);
            
            %prepare for next location (for future code)
            count_eq = count_eq + len_eq_curr;
            count_ineq = count_ineq + len_ineq_curr;                                             
        end
                    
  

    end
end

