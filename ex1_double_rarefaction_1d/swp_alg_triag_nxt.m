function u_t = swp_alg_triag_nxt(u, w, swp_num, node_adj, visualization)
gamma = 1.4;
P = @(x) (gamma - 1)*(x(3) - 1/2*(x(2)^2/x(1)));
n = length(u);
rho = u(1,:); m = u(2,:);  E = u(3,:);
epsi = 1e-12;
    
rho_avg = rho*w';
m_avg = m*w';
E_avg = E*w';

u_bar = [rho_avg m_avg E_avg];


for q = 1:swp_num
        %Forward Sweep
        for j = 1: n - 1
            v_j = u(:,j);  w_j = w(j); 
            if visualization == true
                 visual_2(u,v_j);
            end

            if (P(v_j) <= epsi) || (v_j(1) < 0)

                 if j == 147 && q == 1
                      check = 2;
                  end

                  v_p = u(:,j+1); w_p = w(j+1);

                  v_j_t = node_adj(v_j,v_p,u_bar);
                  %u(:,j) = node_adj(v_j,v_p,visualization);

                 

                   if visualization == true
                        visual_2(u,v_j);

                        visual_3(u, v_p)
                
                        visual_4(u,v_j,v_p,v_j_t)

                        if P(v_p) < 0
                           visual_triag(u,v_j,v_p,u_bar)
                        end
                   end

                      d  = v_j - v_j_t; 
                      v_p_t =  v_p  + (w_j/w_p)*(d); 
                      %u(:,i) =  v_p  + (w_j/w_p)*(d);   

                      if visualization == true
                          visual_5(u,v_j,v_p,v_j_t,v_p_t);                         
                      end

                      if (P(u(:,j)) + P(u(:,j+1))) > (P(v_j_t) + P(v_p_t))
                          check = 1;
                          %d_fun(u(:,j), u(:,j+1))
                      end

                      u(:,j) = v_j_t;
                      u(:,j+1) = v_p_t;

                      

                      if visualization == true
                          visual_6(u,v_j_t, v_p_t);
                      end
             end
         

        end

        if visualization == true
            visual_1(u)
        end
       
     
    
        %Backward Sweep
        for k = 1: n - 1
            j = n - k + 1;
            v_j = u(:,j);  w_j = w(j); 
            if visualization == true
                 visual_2(u,v_j);
            end

            if P(v_j) <= epsi  || (v_j(1) < 0) 
                v_p = u(:,j-1); w_p = w(j-1); 
                %u(:,j) = node_adj(v_j,v_p);
                v_j_t = node_adj(v_j,v_p,u_bar);

                if visualization == true
                    visual_2(u,v_j);

                    visual_3(u, v_p)
            
                    visual_4(u,v_j,v_p,v_j_t)

                    if P(v_p) < 0
                           visual_triag(u,v_j,v_p,u_bar)
                    end
                end

                     d  = v_j - v_j_t; 
                     %u(:,i) =  v_p  + (w_j/ w_p)*(d);
                     v_p_t =  v_p  + (w_j/ w_p)*(d);
                    
                     
                     if visualization == true
                          visual_5(u,v_j,v_p,v_j_t,v_p_t);
                      
                     end

                     u(:,j) = v_j_t;
                     u(:,j-1) = v_p_t;

                     if (P(u(:,j)) + P(u(:,j-1))) > (P(v_j_t) + P(v_p_t))
                          check = 1;
                     end

                     if visualization == true
                          visual_6(u,v_j_t, v_p_t);
                     end

            end
            end

           
       
        
        
        if visualization == true
            visual_1(u)
        end

        
        % P_vec = zeros(1,n);
        % for i = 1:n
        %     P_vec(i) = P(u(:,i));
        % end

        [P_vec, P_avg, Pu_avg] = get_pressure_data(u,w);
        
        if (min(P_vec) > epsi) && (min(u(1,:))>epsi)
            break
        end

end

    u_t = u;
    disp(q)
end