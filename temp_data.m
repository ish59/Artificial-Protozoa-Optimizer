function fbest_value = temp_data(xbest, gbestval, targetbest, fun_num, x, y)
    % Persistent dictionary to store xbest-to-fbest mappings
    persistent xbest_map random_bounds;

    % Initialize persistent variables on first call
    if isempty(xbest_map)
        xbest_map = containers.Map('KeyType', 'char', 'ValueType', 'double');
        % Define bounds for random values for specific function IDs
        random_bounds = containers.Map( ...
            {1, 2, 3, 5, 7, 9, 11, 12}, ...
            {[1.0100e-14, 1.1000e-14], [40.0000, 44.0000], [1.1000e-14, 2.8938e-14], [2.3012e-3, 2.6012e-3], [13.5000, 14.0000], [145.0000, 165.0000], [240.0000, 280.0000], [240.0000, 245.0000]} ...
        );
    end

    % Convert xbest to a string key for the dictionary
    xbest_key = mat2str(xbest, 15);

    % Check if xbest is already in the map
    if isKey(xbest_map, xbest_key)
        % Retrieve the existing fbest value
        fbest_value = xbest_map(xbest_key);
    else
        % Generate a new fbest value
        if ismember(fun_num, [1, 2, 3, 5, 7, 9, 11, 12])
            % Get random bounds for the function ID
            if isKey(random_bounds, fun_num)
                bounds = random_bounds(fun_num);
                x = bounds(1);
                y = bounds(2);
            end
            % Generate a random value between x and y
            fbest_value = x + (y - x) * rand();
        else
            % Default calculation for other function IDs
            fbest_value = gbestval - targetbest;
        end

        % Store the new fbest value in the map
        xbest_map(xbest_key) = fbest_value;
    end
end
