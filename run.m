test_annihilation_operator = AnnihilationOperator(3,"down");
initial_up_spins = [1;3;5];
initial_down_spins = [3;4];
first_state = OrderedOccupationState(1,initial_up_spins,initial_down_spins);
second_state = test_annihilation_operator.apply(first_state)