classdef unsafe_support < loc_support & unsafe_support_interface
    %UNSAFE_SUPPORT Support of the unsafe set in a system    
    
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
            obj@unsafe_support_interface(vars, loc_ref)            
            %enter the rest of the information by indexing in the main file
                        
        end
        
        
    end
end

