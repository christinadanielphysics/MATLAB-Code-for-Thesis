classdef OrderedOccupationState
    properties
        Coefficient
        Up_Spins 
        Down_Spins
    end

    methods
        function obj = OrderedOccupationState(Coefficient,Up_Spins,Down_Spins)
            obj.Coefficient = Coefficient;
            obj.Up_Spins = Up_Spins;
            obj.Down_Spins = Down_Spins;
        end
    end
end