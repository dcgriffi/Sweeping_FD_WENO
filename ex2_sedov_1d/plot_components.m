function plot_components(u, roRef, x, xx, figure_num, Method, tFinal, exact)


if exact == true
    figure(figure_num)
    subplot(3,1,1)
    plot(x, u(1,:),'.-')
     hold on
     plot(xx, roRef, 'ro')
    % legend('Approximate','Exact')
    title(num2str(Method) + " for ro at Time " + tFinal)
    legend('Approximate','Exact')
    subplot(3,1,2)
    plot(x, u(2,:),'.-')
    % hold on
    % plot(x, uRef(:,2), 'ro')
    title(num2str(Method) + " for m at Time " + tFinal)
    % legend('Approximate','Exact')
    subplot(3,1,3)
    plot(x, u(3,:),'.-')
    % hold on
    % plot(x, uRef(:,3), 'ro')
    title(num2str(Method) + "for E at Time " + tFinal)
    % legend('Approximate','Exact')
end

if exact == false

    ro = u(1,:);
    m = u(2,:);
    E = u(3,:);

    v = m./ro;
    p = E*(1.4 - 1) - 1/2*ro.*v.^2*(1.4 - 1) ;
    

    figure(figure_num)
    subplot(5,1,1)
    plot(x, u(1,:),'.-')
    title(num2str(Method) + " for rho at Time " + tFinal)
   % legend('Approximate')
    subplot(5,1,2)
    plot(x, u(2,:),'.-')
    title(num2str(Method) + " for m at Time " + tFinal)
   % legend('Approximate')
    subplot(5,1,3)
    plot(x, u(3,:),'.-')
    title(num2str(Method) + "for E at Time " + tFinal)
   % legend('Approximate')
    subplot(5,1,4)
    plot(x, v,'.-')
    title(num2str(Method) + "for v at Time " + tFinal)
   % legend('Approximate')
    subplot(5,1,5)
    plot(x, p,'.-')
    title(num2str(Method) + "for p at Time " + tFinal)
 %   legend('Approximate')
elseif exact == 6

    ro = u(1,:);
    m = u(2,:);
    E = u(3,:);

    v = m./ro;
    p = E*(1.4 - 1) - 1/2*ro.*v.^2*(1.4 - 1) ;
    

    figure(figure_num)
    subplot(5,1,1)
    plot(x, u(1,:),'.-')
 
    hold on;
     a=load('density.dat'); 
    xx=a(:,1); roRef=a(:,2); 
    x_new = [flip(-xx);xx];
    roRef_new = [flip(roRef);roRef];
    plot(x_new,roRef_new);
    title(num2str(Method) + " for rho at Time " + tFinal)
   % legend('Approximate')
    subplot(5,1,2)
    plot(x, u(2,:),'.-')
    title(num2str(Method) + " for m at Time " + tFinal)
   % legend('Approximate')
    subplot(5,1,3)
    plot(x, u(3,:),'.-')
    title(num2str(Method) + "for E at Time " + tFinal)
   % legend('Approximate')
    subplot(5,1,4)
    plot(x, v,'.-')
    title(num2str(Method) + "for v at Time " + tFinal)
   % legend('Approximate')
    subplot(5,1,5)
    plot(x, p,'.-')
    title(num2str(Method) + "for p at Time " + tFinal)
 %   legend('Approximate')
end

end