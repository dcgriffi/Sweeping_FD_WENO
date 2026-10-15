clear all
close all

%Loading Fortran run results
rho =load('den.dat');  rho = rho(:,2)';
v = load('vel.dat');  v = v(:,2)';
p = load('pre.dat');  p = p(:,2)';


load_Ref =load('density.dat');
rhoRef = load_Ref(:,2)';
xx_Ref = load_Ref(:,1)';


%Preparing function inputs
gamma = 1.4; n = length(rho); xL = -2; xR = 2;  dx = (xR - xL)/n;
m = v.*rho;
E = p./(gamma - 1) + 1/2*rho.*v.^2;
u = [rho;m;E];

%uRef = NaN;
xx = NaN(n + 6,1);
for i = 1:n + 6
    index = i - 3;
    xx(i) = xL + index*dx - dx/2;
end
xx = xx(4:n+3);

% n_ex = length(rhoRef)-6;
% xx_Ref = NaN(n_ex + 6,1);
% for i = 1:n_ex + 6
%     index = i - 3;
%     xx_Ref(i) = xL + index*dx - dx/2;
% end

figure_num = 1;
Method = "Fortran FD WENO RK3";
tFinal = 0.001; 
exact = false;

%Plotting components
plot_components(u, rhoRef, xx, xx_Ref, figure_num, Method, tFinal, exact)