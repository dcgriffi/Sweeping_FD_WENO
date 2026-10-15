clear all
close all

problem = 2; % 

pth = pwd;  % figures are saved to ./results_2D-vortex-accuracy

%Loading Fortran run results
 gamma = 1.4; gamma1 = (gamma - 1); t_final = 0.01; use_lim = 2;
rho_load =load('den.dat');  
xx = rho_load(:,1)'; yy = rho_load(:,2)'; 
rho = rho_load(:,3)'; 

n = sqrt(length(rho));
x = reshape(xx,n,n); y = reshape(yy,n,n);
x  = x(:,1); y = y(1,:);

rho = reshape(rho, n,n);
%v = load('vel1.dat');  v = v(:,3)'; v = reshape(v,n,n);
%w = load('vel2.dat');  w = w(:,3)'; w = reshape(w,n,n);

m1_res = load('m1.dat');  m1_res = m1_res(:,3)'; m1_res = reshape(m1_res,n,n);
m2_res = load('m2.dat');  m2_res = m2_res(:,3)'; m2_res = reshape(m2_res,n,n);

p = load('pre.dat');   p = p(:,3)'; p = reshape(p,n,n);

if problem == 1
    rho_exact = @(x,y,t) 1.+0.99*sin(2*pi*(x + y - 2*t));
    m1 = @(x,y,t) rho_exact(x,y,t);
    m2 = @(x,y,t) rho_exact(x,y,t);
    E = @(x,y,t)  1./gamma1 +rho_exact(x,y,t);

elseif problem == 2
    r= @(x,y,t) sqrt((x-5.-t)^2+(y-5.-t)^2);
    !mag=10.05;
    mag = 10.0828;
    rho_exact = @(x,y,t) (1-gamma1*mag^2/8./gamma/pi^2*exp(1.-r(x,y,t)^2))^(1./gamma1);
    m1 = @(x,y,t) rho_exact(x,y,t)*(1.-mag/2./pi*exp((1.-r(x,y,t)^2)/2.)*(y-5.));
    m2 = @(x,y,t) rho_exact(x,y,t)*(1.+mag/2./pi*exp((1.-r(x,y,t)^2)/2.)*(x-5.));
    E = @(x,y,t) rho_exact(x,y,t)^gamma/gamma1+0.5*m1(x,y,t)^2/rho_exact(x,y,t)+0.5*m2(x,y,t)^2/rho_exact(x,y,t) ;
end

v_exact = @(x,y,t) m1(x,y,t)./(rho_exact(x,y,t));  
w_exact = @(x,y,t) m2(x,y,t)./(rho_exact(x,y,t));   
p_exact = @(x,y,t) gamma1*(E(x,y,t) - 0.5*(m1(x,y,t).^2 + m2(x,y,t).^2)./rho_exact(x,y,t)) ;

rho_Ref = NaN(size(rho)); m1_Ref = NaN(size(rho));
m2_Ref = NaN(size(rho)); p_Ref = NaN(size(rho));
for j = 1:n
    for i = 1:n
        rho_Ref(i,j) = rho_exact(x(i), y(j),t_final);
        m1_Ref(i,j)   = m1(x(i), y(j),t_final);
        m2_Ref(i,j)   = m2(x(i), y(j),t_final);
        p_Ref(i,j)   = p_exact(x(i), y(j),t_final);
    end
end
        
Method = "FD WENO"; exact = true;



%Plotting components
%plot_components(use_lim, rho, v, w, p, rho_Ref, v_Ref,w_Ref, p_Ref, x, y, Method, t_final, exact,pth) 

plot_components(use_lim, rho, m1_res, m2_res, p, rho_Ref, m1_Ref,m2_Ref, p_Ref, x, y, Method, t_final, exact,pth) 

n  = length(rho);
n2 = floor(n/2);

den_approx_slice = [y(:), rho(:, n2)];
den_exact_slice  = [y(:), rho_Ref(:, n2)];
pres_approx_slice = [y(:), p(:, n2)];
pres_exact_slice  = [y(:), p_Ref(:, n2)];

%Save each variable to .dat files
save('den_approx_slice.dat', 'den_approx_slice', '-ascii');
save('den_exact_slice.dat',  'den_exact_slice',  '-ascii');
save('pres_approx_slice.dat', 'pres_approx_slice', '-ascii');
save('pres_exact_slice.dat',  'pres_exact_slice',  '-ascii');