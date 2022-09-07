classdef AnnihilationOperator
    properties
        Spatial_Orbital_Index
        Spin
    end

    methods
        function obj = AnnihilationOperator(Spatial_Orbital_Index,Spin)
            obj.Spatial_Orbital_Index = Spatial_Orbital_Index;
            obj.Spin = Spin;
        end

        function [coefficient,spins] = inspect(obj,initial_state)
            
            coefficient = 0; % annihilation

            initial_spins = [];
            spins = [];
            skips = 0;
            end_length = 0;

            if obj.Spin == "up"
                end_length = length(initial_state.Up_Spins) - 1;
                initial_spins = initial_state.Up_Spins;
                spins = zeros(end_length,1);
            else
                skips = length(initial_state.Up_Spins);
                end_length = length(initial_state.Down_Spins) - 1;
                initial_spins = initial_state.Down_Spins;
                spins = zeros(end_length,1);
            end
            
            counter = 1;
            for extra_skips = 1:length(initial_spins)
                if initial_spins(extra_skips) == obj.Spatial_Orbital_Index
                    skips = skips + extra_skips - 1;
                    coefficient = initial_state.Coefficient * (-1)^skips;
                else
                    if counter ~= (end_length + 1)
                        spins(counter) = initial_spins(extra_skips);
                        counter = counter + 1;
                    end
                end
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