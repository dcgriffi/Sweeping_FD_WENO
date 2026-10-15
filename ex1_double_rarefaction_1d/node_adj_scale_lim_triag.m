function u_1 = node_adj_scale_lim_triag(v_1,v_2,v_bar)
    gamma = 1.4;
    P = @(x) (gamma - 1)*(x(3) - 1/2*(x(2)^2/x(1)));
    d = @(v_1,v_2) sqrt((v_2(3) - v_1(3))^2 + (v_2(2) - v_1(2))^2 + (v_2(1) - v_1(1))^2);
    S = @(t) (1 - t)*v_1 + t*v_2;
    %eps = 1e-1;
   %delta = 1e-6;
    epsi = 1e-13;

    if (P(v_2) > epsi)
        t = (P(v_1)-epsi)/(P(v_1) - P(v_2));
        sprintf('positive neighbor')
    else
        %angle = get_angle(v_1,v_bar',v_2,d);
        t_1 = (P(v_1)-epsi)/(P(v_1) - P(v_bar'));
        % if angle < 180 - eps
        t = min(t_1*d(v_1,v_bar')/d(v_1,v_2), 0.25);

        sprintf('negative pair')

        %t = min(t_1*d(v_1,v_bar')/d(v_1,v_2),3/4);

        % else
        % t = t_1;1/2
        % S = @(t) (1 - t)*v_1 + t*v_bar';
        % end
    end



    u_1 = S(t);


    %Need to handle case where t < 0
    % if t <= 0 || t > 1  
    %     u_1 = v_1;
    % end
 
   end