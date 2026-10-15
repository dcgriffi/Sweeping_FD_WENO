function plot_components(use_lim, rho, v, w, p, rho_Ref, v_Ref,w_Ref, p_Ref, x, y, Method, t_final, exact,pth)

n = length(rho);
n2 = floor(n/2);
if exact == true

    figure(1)
    mesh(x,y,rho)
    title(num2str(Method) + " for rho Approx at Time " + t_final+ " with Lim " + use_lim+ " and n = " + n)
    view(20,20)

    saveas(1, sprintf('%s/results_2D-vortex-accuracy/euler_v_approx_lim_%d_tFinal_%d_plot.jpg',pth,use_lim,t_final),'jpeg')
    
    figure(2)
    mesh(x,y,rho_Ref)
    title(num2str(Method) + "  for rho Exact at Time " + t_final+ " with Lim " + use_lim+ " and n = " + n)
    view(20,20)

    %% Density Plot
figure
plot(y, rho(:,n2),'.-','LineWidth',1.2, 'MarkerSize',8)
hold on
plot(y, rho_Ref(:,n2), 'ro','MarkerSize',6,'LineWidth',1.2)
title("Density","FontSize",14)
xlabel("y","FontSize",12)
ylabel("Density","FontSize",12)
legend('Approximate','Exact','Location','southoutside','Orientation','horizontal','FontSize',10)
grid on
set(gca,'FontSize',12,'LineWidth',1.1)
%saveas(gcf, sprintf('%s/results_2D-vortex-accuracy/density_lim_%d_tFinal_%d.jpg', ...
 %     pth, use_lim, t_final), 'jpeg')


%% Pressure Plot
figure
plot(y, p(:,n2),'.-','LineWidth',1.2,'MarkerSize',8)
hold on
plot(y, p_Ref(:,n2), 'ro','MarkerSize',6,'LineWidth',1.2)
title("Pressure","FontSize",14)
xlabel("y","FontSize",12)
ylabel("Pressure","FontSize",12)
legend('Approximate','Exact','Location','southoutside','Orientation','horizontal','FontSize',10)
grid on
set(gca,'FontSize',12,'LineWidth',1.1)
%saveas(gcf, sprintf('%s/results_2D-vortex-accuracy/pressure_lim_%d_tFinal_%d.jpg', ...
 %     pth, use_lim, t_final), 'jpeg')

% figure(3)
% subplot(2,1,1)
% plot(y, rho(:,n2),'.-')
% hold on
% plot(y, rho_Ref(:,n2), 'ro')
% title("Density")
% legend('Approximate','Exact','Location','southoutside','Orientation','horizontal')
% 
% subplot(2,1,2)
% plot(y, p(:,n2),'.-')
% hold on
% plot(y, p_Ref(:,n2), 'ro')
% title("Pressure")
% legend('Approximate','Exact','Location','southoutside','Orientation','horizontal')
% 
% saveas(3, sprintf('%s/results_2D-vortex-accuracy/euler_slices_lim_%d_tFinal_%d_plot_x.jpg', ...
%       pth, use_lim, t_final), 'jpeg')


% figure(4)
% subplot(2,1,1)
% plot(x, squeeze(rho(n2,:)),'.-')
% hold on
% plot(x, squeeze(rho_Ref(n2,:)), 'ro')
% title("Density")
% legend('Approximate','Exact','Location','southoutside','Orientation','horizontal')
% 
% subplot(2,1,2)
% plot(x, squeeze(p(n2,:)),'.-')
% hold on
% plot(x, squeeze(p_Ref(n2,:)), 'ro')
% title("Pressure")
% legend('Approximate','Exact','Location','southoutside','Orientation','horizontal')

    %saveas(2, sprintf('%s/results_2D-vortex-accuracy/euler_v_exact_lim_%d_tFinal_%d_plot.jpg',pth,use_lim,t_final),'jpeg')

    % figure(3)
    % subplot(4,1,1)
    % plot(y, rho(:,n2),'.-')
    % hold on
    % plot(y, rho_Ref(:,n2), 'ro')
    % title(num2str(Method) + " at x = (xR - xL)/2 for den at Time " + t_final+ " with Lim " + use_lim+ " and n = " + n)
    % legend('Approximate','Exact')
    % subplot(4,1,2)
    % plot(y, v(:,n2),'.-')
    % hold on
    % plot(y, v_Ref(:,n2), 'ro')
    % title(num2str(Method) + " at x = (xR - xL)/2 for m1 at Time " + t_final+ " with Lim " + use_lim+ " and n = " + n)
    % legend('Approximate','Exact')
    % subplot(4,1,3)
    % plot(y, w(:,n2),'.-')
    % hold on
    % plot(y, w_Ref(:,n2), 'ro')
    % title(num2str(Method) + " at x = (xR - xL)/2 for m2 at Time " + t_final+ " with Lim " + use_lim+ " and n = " + n)
    % legend('Approximate','Exact')
    % subplot(4,1,4)
    % plot(y, p(:,n2),'.-')
    % hold on
    % plot(y, p_Ref(:,n2), 'ro')
    % title(num2str(Method) + " at x = (xR - xL)/2 for pres at Time " + t_final + " with Lim " + use_lim+ " and n = " + n)
    % legend('Approximate','Exact')
    % 
    % saveas(3, sprintf('%s/results_2D-vortex-accuracy/euler_slices_lim_%d_tFinal_%d_plot_x.jpg',pth,use_lim,t_final),'jpeg')
    % 
    % 
    % figure(4)
    % subplot(4,1,1)
    % plot(x, squeeze(rho(n2,:)),'.-')
    % hold on
    % plot(x, squeeze(rho_Ref(n2,:)), 'ro')
    % title(num2str(Method) + " at y = (yR - yL)/2 for den at Time " + t_final+ " with Lim " + use_lim + " and n = " + n)
    % legend('Approximate','Exact')
    % subplot(4,1,2)
    % plot(x, squeeze(v(n2,:)),'.-')
    % hold on
    % plot(x, squeeze(v_Ref(n2,:)), 'ro')
    % title(num2str(Method) + " at y = (yR - yL)/2 for m1 at Time " + t_final+ " with Lim " + use_lim+ " and n = " + n)
    % legend('Approximate','Exact')
    % subplot(4,1,3)
    % plot(x, squeeze(w(n2,:)),'.-')
    % hold on
    % plot(x, squeeze(w_Ref(n2,:)), 'ro')
    % title(num2str(Method) + " at y = (yR - yL)/2 for m2 at Time " + t_final+ " with Lim " + use_lim+ " and n = " + n)
    % legend('Approximate','Exact')
    % subplot(4,1,4)
    % plot(x, squeeze(p(n2,:)),'.-')
    % hold on
    % plot(x, squeeze(p_Ref(n2,:)), 'ro')
    % title(num2str(Method) + " at y = (yR - yL)/2 for pres at Time " + t_final+ " with Lim " + use_lim+ " and n = " + n)
    % legend('Approximate','Exact')

   % saveas(4, sprintf('%s/results_2D-vortex-accuracy/euler_slices_lim_%d_tFinal_%d_plot_y.jpg',pth,use_lim,t_final),'jpeg')


end

if exact == false

    %Error and convergence analysis
    %Plotting at desired time
    figure(3)
    subplot(4,1,1)
    plot(x, diag(rho),'.-')
    hold on
    title(num2str(Method) + " at x=y for ro at Time " + t_final)
    legend('Approximate')
    subplot(4,1,2)
    plot(x, diag(v),'.-')
    hold on
    title(num2str(Method) + " at x=y for vel1 at Time " + t_final)
    legend('Approximate')
    subplot(4,1,3)
    plot(x, diag(w),'.-')
    hold on
    title(num2str(Method) + "with x=y for vel2 at Time " + t_final)
    legend('Approximate')
    subplot(4,1,4)
    plot(x, diag(p),'.-')
    hold on
    title(num2str(Method) + "with x=y for pres at Time " + t_final)
    legend('Approximate')


end