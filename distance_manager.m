classdef distance_manager < manager_interface
    %DISTANCE_MANAGER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
 
        supp;
        meas_wass;
        
        dist;
    end
    
    methods
        function obj = distance_manager(unsafe_supp, f)
            %DISTANCE_MANAGER Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;

            %set up the location
            loc_curr = location(1, unsafe_supp, f, []);
            obj@manager_interface(loc_curr);

            %information about unsafe set and distance
            obj.supp = unsafe_supp;
            
            %create wasserstein measure
            %assume a single unsafe set, 
            %TODO: generalize to multiple later
            supp_joint = [obj.supp.X; obj.supp.X_unsafe];
            
            obj.meas_wass = meas_joint();
                        
        end
        
  %% Formulating and solving program
        
        function [objective, mom_con, supp_con] = cons(obj,d, Tmax)
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
            [~, cons_eq, cons_ineq] = obj.loc.all_cons(d);
            
            %TODO: joint marginal constraint
            %wasserstein objective
            
            %finalize moment constraints
            
            %mass of initial measure sums to one
            mass_init_con = (obj.loc.mass_init() - 1 == 0);
            
            if obj.loc.supp.TIME_INDEP
                mass_occ_con = (obj.loc.mass_occ() <= Tmax);
            else
                mass_occ_con = [];
            end

            mom_con = [cons_eq; cons_ineq; mass_init_con; mass_occ_con];
            

        end
    end
end

