function [P_vec, P_avg, Pu_avg] = get_pressure_data(u,weights)
        P = @(x) x(3) - 1/2*(x(2)^2/x(1));
        n = length(u);
        
        P_vec = zeros(1,n);
        for i = 1:n
            P_vec(i) = P(u(:,i));
        end

        rho = u(1,:); m = u(2,:);  E = u(3,:);
        
        rho_avg = rho*weights';
        m_avg = m*weights';
        E_avg = E*weights';
        
        u_avg = [rho_avg m_avg E_avg];
     
        P_avg = P_vec*weights';
        Pu_avg = P(u_avg);
end
