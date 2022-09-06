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

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end