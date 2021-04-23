classdef unsafe_support < loc_support
    %UNSAFE_SUPPORT Support of the unsafe set in a system    
    
    properties
%         vars = struct('x', [], 'y', []);        
        
        %unsafe set (y)
        
        X_unsafe = [];
        
        dist; 
    end
    
    methods
        function obj = unsafe_support(vars, loc_ref)
            %UNSAFE_SUPPORT Construct an instance of this class
            %   Detailed explanation goes here
%             obj.Property1 = inputArg1 + inputArg2;
        

            %TODO: redo interface with 'distance' function
               
            if nargin == 1
                loc_ref = [];
            end
            
            obj@loc_support(vars, loc_ref);
            
            if nargin > 1                
                obj.X_unsafe = subs_vars(loc_ref.X_unsafe, loc_ref.vars.y, obj.vars.y);                
                obj.dist = loc_ref.dist;                   
            end
            
            %enter the rest of the information by indexing in the main file
                        
        end
        
        function vars_out = get_vars_joint(obj)
            vars_out = [obj.vars.x;
                        obj.vars.y];
        end
        
        function Xu = get_X_unsafe(obj)
            Xu = obj.X_unsafe;
        end
        
        
    end
end
