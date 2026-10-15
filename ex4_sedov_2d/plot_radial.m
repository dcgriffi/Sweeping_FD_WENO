clear all
close all

problem = 2; % 


%Loading Fortran run results
 gamma = 1.4; gamma1 = (gamma - 1); t_final = 1; use_lim = 2;
rho_load =load('den_mat.dat');  
xx = rho_load(:,1)'; yy = rho_load(:,2)'; 
rho = rho_load(:,3)'; 

n = sqrt(length(rho));
x = reshape(xx,n,n); y = reshape(yy,n,n);
x  = x(:,1); y = y(1,:);

rho = reshape(rho, n,n);


Method = "FD WENO"; exact = true;



%Plotting components
%plot_components(use_lim, rho, v, w, p, rho_Ref, v_Ref,w_Ref, p_Ref, x, y, Method, t_final, exact,pth) 

%plot_components(use_lim, rho, m1_res, m2_res, p, rho_Ref, m1_Ref,m2_Ref, p_Ref, x, y, Method, t_final, exact,pth) 

%Load exact solution

figure(1)
mesh(x,y,rho)
title(num2str(Method) + " for rho Approx at Time " + t_final+ " with Lim " + use_lim+ " and n = " + n)
view(20,20)

exact_sol =  load('sedov.dat');
r = exact_sol(:,1);
v_Ref = exact_sol(:,2);
den_Ref = exact_sol(:,3);
p_Ref = exact_sol(:,4);
u_Ref = exact_sol(:,5);

figure(7)
plot(r, den_Ref)
hold on 
plot(x, rho(1,:),'x', 'MarkerSize', 8, 'LineWidth', 1.5);
legend('Exact', 'y = 0')

%Radial projection
X = x; Y = y;
% Center at (0,0) for quarter-domain symmetry
R = sqrt(X.^2 + Y.^2);

% Radial bins
%nBins = 200;
nBins = 640;
edges = linspace(0, 1.3, nBins+1);
rMid  = 0.5 * (edges(1:end-1) + edges(2:end));

% Flatten arrays for binning
R_flat = R(:);
F_flat = rho(:);


% Bin assignment
[~,~,bin] = histcounts(R_flat, edges);

% Remove out-of-range points (bin == 0)
mask = bin > 0;
bin  = bin(mask);
F_flat = F_flat(mask);

% Radial averaging
profile = accumarray(bin, F_flat, [nBins 1], @mean, NaN);

% Plot
figure(8);
plot(rMid, profile, 'x', 'MarkerSize', 8, 'LineWidth', 1.5);
xlabel('r'); ylabel('f(r)');
title('Density');
grid on;
hold on;
plot(r, den_Ref)
legend('Radial Projection', 'Exact')


figure(9)
plot(r, den_Ref)
hold on 
plot(x,rho(:,1),'x', 'MarkerSize', 8, 'LineWidth', 1.5);
legend('Exact', 'y = 0')