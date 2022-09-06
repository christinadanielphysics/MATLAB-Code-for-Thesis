test_creation_operator = CreationOperator(2,"up");
initial_up_spins = [1];
initial_down_spins = [];
first_state = OrderedOccupationState(1,initial_up_spins,initial_down_spins);
second_state = test_creation_operator.apply(first_state)
second_state.Up_Spins;