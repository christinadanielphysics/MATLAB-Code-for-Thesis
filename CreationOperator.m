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

        function [coefficient,spins] = inspect(obj,initial_state)
     
            initial_spins = [];
            spins = [];
            skips = 0;
            operator_insertion_index = 1;
            end_length = 1;

            if obj.Spin == "up"
                operator_insertion_index = length(initial_state.Up_Spins) + 1;
                end_length = length(initial_state.Up_Spins) + 1;
                initial_spins = initial_state.Up_Spins;
                spins = zeros(end_length,1);
            else
                skips = length(initial_state.Up_Spins);
                operator_insertion_index = length(initial_state.Down_Spins) + 1;
                end_length = length(initial_state.Down_Spins) + 1;
                initial_spins = initial_state.Down_Spins;
                spins = zeros(end_length,1);
            end

            for extra_skips = 1:length(initial_spins)
                if initial_spins(extra_skips) == obj.Spatial_Orbital_Index
                    coefficient = 0;
                    spins = [];
                    return
                elseif initial_spins(extra_skips) > obj.Spatial_Orbital_Index
                    operator_insertion_index = extra_skips;
                    break;
                else
                    skips = skips + 1;
                    spins(extra_skips) = initial_spins(extra_skips);
                    operator_insertion_index = extra_skips + 1;
                end
            end

            coefficient = initial_state.Coefficient * (-1)^skips;
            spins(operator_insertion_index) = obj.Spatial_Orbital_Index;
            
            while operator_insertion_index < end_length
                spins(operator_insertion_index+1) = initial_spins(operator_insertion_index);
                operator_insertion_index = operator_insertion_index + 1;
            end

            return
            
        end

        function final_state = apply(obj,initial_state)

            if initial_state.Coefficient == 0
                final_state = OrderedOccupationState(0,[],[]);
                return
            end

            final_up_spins = initial_state.Up_Spins;
            final_down_spins = initial_state.Down_Spins;
            final_coefficient = initial_state.Coefficient;

            if obj.Spin == "up"
                [final_coefficient,final_up_spins] = inspect(obj,initial_state);
            else
                [final_coefficient,final_down_spins] = inspect(obj,initial_state);
            end

            final_state = OrderedOccupationState(final_coefficient,final_up_spins,final_down_spins);
            return
        end
    end
end

