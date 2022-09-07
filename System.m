classdef System

    properties
        Number_of_Spatial_Orbitals
        Number_of_Spin_Up_Electrons 
        Number_of_Spin_Down_Electrons
    end

    methods
        function obj = System(Number_of_Spatial_Orbitals,Number_of_Spin_Up_Electrons,Number_of_Spin_Down_Electrons)
            obj.Number_of_Spatial_Orbitals = Number_of_Spatial_Orbitals;
            obj.Number_of_Spin_Up_Electrons = Number_of_Spin_Up_Electrons;
            obj.Number_of_Spin_Down_Electrons = Number_of_Spin_Down_Electrons;
        end

        function basis_states = get_basis_states(obj)
            up_array_matrix = nchoosek(1:obj.Number_of_Spatial_Orbitals,obj.Number_of_Spin_Up_Electrons);
            down_array_matrix = nchoosek(1:obj.Number_of_Spatial_Orbitals,obj.Number_of_Spin_Down_Electrons);

            [numRows_for_up,numCols_for_up] = size(up_array_matrix);
            [numRows_for_down,numCols_for_down] = size(down_array_matrix);

            basis_states = OrderedOccupationState.empty;

            counter = 1;

            for up_row = 1:numRows_for_up

                up_array = zeros(numCols_for_up,1);
                for up_col = 1:numCols_for_up
                    up_array(up_col,1) = up_array_matrix(up_row,up_col);
                end

                for down_row = 1:numRows_for_down

                    down_array = zeros(numCols_for_down,1);
                    for down_col = 1:numCols_for_down
                        down_array(down_col,1) = down_array_matrix(down_row,down_col);
                    end

                    state = OrderedOccupationState(1,up_array,down_array);
                    basis_states(1,counter) = state;
                    counter = counter + 1;
                end
            end

            return

        end

    end
end