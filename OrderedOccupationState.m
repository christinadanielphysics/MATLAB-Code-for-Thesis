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
                up_string = up_string + string(obj.Up_Spins(element)) + "_";
            end

            down_string = "";
            for element = 1:length(obj.Down_Spins)
                down_string = down_string + string(obj.Down_Spins(element)) + "_";
            end

            return
        end

        function product = scalar_product(obj,right_state)
            product = 0;
            if isequal(obj.Up_Spins,right_state.Up_Spins) && isequal(obj.Down_Spins,right_state.Down_Spins)
                product = obj.Coefficient * right_state.Coefficient;
            end
            return
        end

        function number = scalar_product_one_with_many(obj,final_states)
            number = 0;
            for index = 1:length(final_states)
                if final_states(index).Coefficient ~= 0
                    number = number + obj.scalar_product(final_states(index));
                end
            end
            return
        end

    end
end