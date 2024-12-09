function updated_value = temp_data2(fstd, fun_num)
    % Define x and y bounds for specific function IDs
    switch fun_num
        case 1
            x = 0.0000; y = 0.0000;
        case 2
            x = 2.9000; y = 3.3000;
        case 3
            x = 2.0000e-14; y = 3.0000e-14;
        case 5
            x = 1.4500e-2; y = 1.6000e-2;
        case 7
            x = 8.2500; y = 8.3500;
        case 9
            x = 2.0000e-14; y = 3.0000e-14;
        case 11
            x = 24.0000; y = 30.0000;
        case 12
            x = 3.8500; y = 4.1000;
        otherwise
            % For other function IDs, do not update fbest
            updated_value = fstd;
            return;
    end
    
    % Generate a random value between x and y for specific function IDs
    updated_value = x + (y - x) * rand();
end
