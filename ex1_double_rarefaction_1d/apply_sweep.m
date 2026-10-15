function uHold = apply_sweep(uHold,weights, swp_num, swp_alg, node_adj, visualization)
    [P_vec, P_avg, Pu_avg] = get_pressure_data(uHold,weights);
                         
    if min(uHold(1,:)) <= 1e-13
        uHold(1,:) = scalar_pos(uHold(1,:), weights);
    end

    [P_vec, P_avg, Pu_avg] = get_pressure_data(uHold,weights);


    if Pu_avg < 0
        check = 3;
    end

    if (min(P_vec)<=1e-13 && Pu_avg > 0) 
        uHold(1,:) = scalar_pos(uHold(1,:), weights);
        [P_vec, P_avg, Pu_avg] = get_pressure_data(uHold,weights);

        if Pu_avg > 1e-13
            u_t = swp_alg(uHold, weights, swp_num, node_adj, visualization);
            uHold = u_t;
        else
            check = 2;
        end
    end
end