classdef CreationOperator
    properties
        Spatial_Orbital_Index
        Spin
    end

    methods
        function obj = CreationOperator(Spatial_Orbital_Index,Spin)
            obj.Spatial_Orbital_Index = Spatial_Orbital_Index;
            obj.Spin = Spin;
        end

        function final_state = apply(obj,initial_state)
            import OrderedOccupationState
            if initial_state.Coefficient == 0
                final_state = OrderedOccupationState(0,int16.empty(0,0),int16.empty(0,0));
                return
            end
        end
    end
end

