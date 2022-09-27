classdef CompressiveSensing

    properties
        n
        y
        perm
        zero
        compressive_threshold
    end

    methods
        function obj = CompressiveSensing(n,y,perm,zero,compressive_threshold)
            obj.n = n;
            obj.y = y;
            obj.perm = perm;
            obj.zero = zero;
            obj.compressive_threshold = compressive_threshold;
        end

        function s = compute(obj)
            % Theta*s == y from Steve Brunton is the same as A*x = b from Emmanuel Candes
            n = obj.n;
            y = obj.y;

            Psi = dct(eye(n,n)); % build Psi
            Theta = Psi(obj.perm,:); % measure rows of Psi
            
            % L1-Minimization using CVX
            cvx_begin;
                variable s(n);
                minimize( norm(s,1) );
                subject to
                    Theta*s == y;
            cvx_end;
        end

        function [combined_w_values,combined_weights] = combine(obj,input_w_values,input_weights)
            number_of_points = length(input_weights);
            combined_w_values = [];
            combined_weights = [];
            counter = 1;
            while counter <= number_of_points
                initial_counter = counter;
                weight = input_weights(counter);
                w_value = input_w_values(counter);
                sum = 0;
                while abs(input_weights(counter)) > obj.zero && counter ~= number_of_points
                    sum = sum + input_weights(counter); % add the weights
                    weight = sum; % update the weight
                    w_value = 0.5 * (w_value + input_w_values(counter)); % average the angular frequencies
                    counter = counter + 1;
                end
                combined_w_values(initial_counter) = w_value;
                combined_weights(initial_counter) = weight;
                counter = counter + 1;
            end
        end

        function [new_w_values,new_weights] = chop(obj,combined_w_values,combined_weights)
            array_length = length(combined_weights);
            new_w_values = [];
            new_weights = [];
            counter = 1;
            for index = 1:array_length
                if abs(combined_weights(index))/max(abs(combined_weights)) >= obj.compressive_threshold
                    new_w_values(counter) = combined_w_values(index);
                    new_weights(counter) = combined_weights(index);
                    counter = counter + 1;
                end
            end
        end
    end
end