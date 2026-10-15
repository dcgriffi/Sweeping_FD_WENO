function new_arr = scalar_pos(arr, weights)
    % BND_PRS_SWP - Apply convservative bound preserving sweep to input array.
    % 
    %   Arguments:
    %      arr - is an array which indicates number of weights
    %      weights - is an array with positive weights summing to one and
    %         same dimensions as arr
    %      m - is a desired lower bount for the array
    %      M - is a a desired upper bound for the array
    %
    %   Returns:
    %      new_arr - is an updated array with values between m and M that
    %      conserves sum of weights with array value

    n = length(arr); m = 1e-16;
    
    % TO DO:
    %Check that array sizes match
    %Check that weights sum to one

    avg_array = sum(arr)/n;

    for j = 1: n - 1
        v_j = arr(j); v_jp1 = arr(j + 1); w_j = weights(j); w_jp1 = weights(j + 1);
        if v_j < m
          arr(j + 1) =  v_jp1  + (w_j/w_jp1)*(v_j - m);
          arr(j) = m;
          if (arr(j) + arr(j+1)) ~= (v_j + v_jp1)
              check = 4;
          end
        end
    end

    for k = 1: n - 1
        %j = n:2
        j = n - k + 1;
        v_j = arr(j); v_jm1 = arr(j - 1); w_j = weights(j); w_jm1 = weights(j - 1);
        if v_j < m
          arr(j - 1) =  v_jm1  + (w_j/ w_jm1)*(v_j - m);
          arr(j) = m;
        end
    end

    new_arr = arr;

    avg_array = sum(arr)/n;

end