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

        function [up_string,down_string] = get_strings(obj)
            up_string = "";
            for element = 1:length(obj.Up_Spins)
                up_string(element) = string(obj.Up_Spins(element)) + "_";
            end

            down_string = "";
            for element = 1:length(obj.Down_Spins)
                down_string(element) = string(obj.Down_Spins(element)) + "_";
            end

            return
        end

    end
end