function [ro,m,E,v,p] = get_variables(u, gam)

ro = u(1,:);
m = u(2,:);
E = u(3,:); 

v = m./ro;
p = (gam - 1)*(E - 1/2*(m.^2)./ro) ;
end


    % %Build matrix R for each point
        % R = zeros(3,3,nSteps+1);
        % for kk = 1:nSteps+1
        %     j = kk + 2;
        % 
        %     t2_j = t2(j); qm_j = qm(j);
        %     Vxm_j = Vxm(j); Hxm_j = Hxm(j); cm_j = cm(j); 
        % 
        %     R(1,:,kk) = 1;
        %     R(2,1,kk) = Vxm_j - cm_j;
        %     R(2,2,kk) = Vxm_j;
        %     R(2,3,kk) = Vxm_j + cm_j;
        %     R(3,1,kk) = Hxm_j - t2_j;
        %     R(3,2,kk) = qm_j;
        %     R(3,3,kk) = Hxm_j + t2_j;
        % end
