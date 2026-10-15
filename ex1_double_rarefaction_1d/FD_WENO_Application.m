clear all;
close all
format short e;

%1D FD WENO
    %Solving the Euler equations with FD WENO
    %using the Lax Friedrichs method.

%Choose Problem Number
problem =5; %Pick 1,2,3,4,5,6

%Fifth Order WENO Scheme
Method = "Fifth Order WENO ";

order = 3; 

sweeping = true;


%Comment or uncomment if you want the shock for T=2.
%avoidShock = true; %(without shock)
avoidShock = false; %(with shock)

fig_num = 0;

if problem == 1 %Linear
    tFinal = 2;

    gamma = 1.4;
    roInitial = @(x) 1+ 0.2*sin(x); 
    vInitial = @(x) 1;
    pInitial = @(x) 1;
    xL = 0; xR = 2*pi;

    periodic = true;
    exact = true;
end

if problem == 2 %Non-linear
    tFinal = 0.3;

    gamma = 3;
    roInitial = @(x) 1 + 1/2*sin(x); 
    vInitial = @(x) 1/2*sin(x);
    pInitial = @(x) (1 + 1/2*sin(x)).^gamma;
    xL = 0; xR = 2*pi;  

    periodic = true;
    exact = true;
end

if problem == 3 %Shock tube
    %tFinal = 2;
    tFinal = 0.4;

    gamma = 1.4;
    idx1 = @(x) x <=0;
    idx2 = @(x) x > 0;
    roInitial = @(x) zeros(1,length(x))+ 1*idx1(x) + 0.125*idx2(x); 
    vInitial = @(x) zeros(1,length(x));
    pInitial = @(x) zeros(1,length(x))+ 1*idx1(x) + 0.1*idx2(x);
    xL = -5; xR = 5;
     xL = -1; xR = 1;

    periodic = false;
    exact = false;
    uRef = NaN;
end

if problem == 4 %Lax
    gamma = 1.4;
    tFinal = 1.3;
    idx1 = @(x) x <=0;
    idx2 = @(x) x > 0;
    roInitial = @(x) zeros(1,length(x))+ 0.445*idx1(x) + 0.5*idx2(x); 
    vInitial = @(x) zeros(1,length(x))+ 0.698*idx1(x) + 0*idx2(x);
    pInitial = @(x) zeros(1,length(x))+ 3.528*idx1(x) + 0.571*idx2(x);
    xL = -5; xR = 5;

    periodic = false;
    exact = false;
    uRef = NaN;
end

if problem == 5 %Double Rarefaction (1,2,3 problem) density 
    gamma = 1.4;
    tFinal =0.3;
    idx1 = @(x) x <=0;
    idx2 = @(x) x > 0;
    roInitial = @(x) zeros(1,length(x))+ 7*idx1(x) + 7*idx2(x); 
    vInitial = @(x) zeros(1,length(x))+ -1*idx1(x) + 1*idx2(x);
    pInitial = @(x) zeros(1,length(x))+ 0.2*idx1(x) + 0.2*idx2(x);
     % roInitial = @(x) zeros(1,length(x))+ 1*idx1(x) + 1*idx2(x); 
     % mInitial = @(x) zeros(1,length(x))+ -1*idx1(x) + 1*idx2(x);
     % vInitial = @(x) mInitial(x)./roInitial(x);
     % EInitial = @(x) zeros(1,length(x))+ 5*idx1(x) + 5*idx2(x);
     % pInitial = @(x) (gamma - 1)*roInitial(x).*EInitial(x);
    xL = -0.5; xR = 0.5;

    cfl = 0.85;
    periodic = false;
    exact = false;
    uRef = NaN;
end


if problem == 6 %1D Sedov Blast
    gamma = 1.4;
    tFinal = 0.0006;
    idx1 = @(x) x >= (-2/402);
    idx2 = @(x) x > (2/402);
    roInitial = @(x) ones(1,length(x)); 
    vInitial = @(x) zeros(1,length(x));
    EInitial = @(x) zeros(1,length(x)) + (3200000*201)*idx1(x) - (3200000*201)*idx2(x) + 1e-12*ones(1,length(x));
    pInitial = @(x) (gamma - 1)*(EInitial(x) - 1/2*roInitial(x).*vInitial(x).^2);
    xL = -3; xR = 3;

    periodic =true;
    exact = false;
    uRef = NaN;
end

if problem == 6 %1D Sedov Blast
    gamma = 1.4;
    tFinal = 0.0004;
    idx1 = @(x) x >= (-2/402);
    idx2 = @(x) x > (2/402);
    roInitial = @(x) ones(1,length(x)); 
    vInitial = @(x) zeros(1,length(x));
    EInitial = @(x) zeros(1,length(x)) + (3200000*201)*idx1(x) - (3200000*201)*idx2(x) + 1e-12*ones(1,length(x));
    pInitial = @(x) (gamma - 1)*(EInitial(x) - 1/2*roInitial(x).*vInitial(x).^2);
    xL = -3; xR = 3;

    periodic =true;
    exact = false;
    uRef = NaN;
end

if problem == 7 %High Mach number astrophysical jet
    gamma = 1.4;
    tFinal = 0.00001;
    %idx1 = @(x) x >= (-1/402);
    %idx2 = @(x) x > (1/402);
    roInitial = @(x) ones(1,length(x))*5; 
    vInitial = @(x) zeros(1,length(x));
    %EInitial = @(x) ones(1,length(x))*5;
    pInitial = @(x) ones(1,length(x))*0.4127;
    xL = -0.05; xR = 0.05;

    periodic =false;
    exact = false;
    uRef = NaN;
end


%Generate points for exact solution (only do once)
x0 = linspace(pi/400,pi-pi/400,3657); 

%Range of points on x for convergence study
if problem == 1 || problem == 2
    nArray=[20,40,80,160,320]; 
elseif problem == 3 || problem == 4 || problem == 5
    nArray = 200;
elseif problem == 6
    nArray = 801;
end

%Preallocate arrays for memory
dx = zeros(length(nArray),1); dt = zeros(length(nArray),1); 
L2Error = zeros(length(nArray),1); L1Error = zeros(length(nArray),1); 
LinfError= zeros(length(nArray),1); err = zeros(length(nArray),1); 

%Cycle through all ranges of dx
for kk=1:length(nArray)

    %Build spatial array
    nSteps = nArray(kk);
    dx(kk)=(xR-xL)/nSteps;
    h=dx(kk); 
    x = (xL:h:xR); 
    x = (x(1:end-1)+x(2:end))/2; %Evaluating all points at xAvg
    
    
    ro = roInitial(x);
    %E = EInitial(x);
    v = vInitial(x);
    p = pInitial(x);

    %plot(x,E)


   
    if problem ==1
        v = v*ones(1,nSteps);
    end

    if problem ==1 
        p = p*ones(1,nSteps);
    end

    m = v.*ro;
    E = p./(gamma - 1) + 1/2*ro.*v.^2;

    u = [ro;m;E];

    weights = ones(1,nSteps)*1/nSteps; 

    fu = [m;(ro.*v.^2 + p);v.*(E + p)] ;

    t = 0;



    %Choose options for sweeping   algorithm   
    swp_alg = @swp_alg_triag_nxt;
    swp_num =200000 ;     
    visualization = false;
    node_adj = @node_adj_scale_lim_triag; 
    %node_adj = @node_adj_ortho;   

    uTrack = u;
    while t < tFinal   
        %Compute eigenvalues
         
        ro = u(1,:);
        m = u(2,:);
        E = u(3,:);
    
        v = m./ro;
        p = E*(gamma - 1) - 1/2*ro.*v.^2*(gamma - 1) ;
        
        fu = [m; (ro.*v.^2 + p); v.*(E + p)];

        c = sqrt((gamma*p)./ro);

        lam1 = v - c;
        lam2 = v;
        lam3 = v + c;

        alpha1 = max(abs(lam1));
        alpha2 = max(abs(lam2));
        alpha3 = max(abs(lam3));
   
        alphas = [alpha1 alpha2 alpha3];
        alpha = max(alphas);

        if order == 3
            k = cfl*(h)/(alpha); 
        elseif order == 5
            k = h^(5/3); %Try to get fifth order acuracy by taking smaller dt
        end

        if tFinal - t < k
             k = tFinal - t;
        end

        l= k/h;
        
        %First update
        [fhatM, fhatP] = FD_WENO_Step(u,gamma,nSteps,periodic);
        uHold = u - l.*(fhatP - fhatM);
        % 
        % uHold(:,1:3) = uTrack(:,1:3);   
        % uHold(:,end-2:end) = uTrack(:,end-2:end);  
        % %uTrack = uHold;

        if sweeping == true
            uHold = apply_sweep(uHold,weights, swp_num, swp_alg, node_adj, visualization);
        end

        %Second update 
        [fhatM, fhatP] = FD_WENO_Step(uHold,gamma,nSteps,periodic);
        uHold = 3/4*u + 1/4*(uHold - l.*(fhatP - fhatM));
       % 
       %  uHold(:,1:3) = uTrack(:,1:3);   
       %  uHold(:,end-2:end) = uTrack(:,end-2:end);  
       % % uTrack = uHold;

        if sweeping == true
            uHold = apply_sweep(uHold,weights, swp_num, swp_alg, node_adj, visualization);
        end


        %Third update
        [fhatM, fhatP] = FD_WENO_Step(uHold,gamma,nSteps,periodic);
        uHold = 1/3*u+ 2/3*(uHold - l.*(fhatP - fhatM));

       %  uHold(:,1:3) = uTrack(:,1:3);   
       %  uHold(:,end-2:end) = uTrack(:,end-2:end);  
       % % uTrack = uHold;

        if sweeping == true
            uHold = apply_sweep(uHold,weights, swp_num, swp_alg, node_adj, visualization);
        end
        
        u = uHold;
        t = t + k;

     
    end

    %Build exact solution
    if problem == 1
        roExact = @(x,t) 1 + 0.2*sin(x - t);
        vExact = @(x,t) 1;
        pExact = @(x,t) 1;

        roRef = roExact(x',tFinal);
        vRef = vExact(x',tFinal);
        pRef = pExact(x',tFinal);
    
        mRef = vRef.*roRef;
        ERef= pRef./(gamma - 1) + 1/2*roRef.*vRef.^2;
    elseif problem ==2
        w = solveBurgers2(-sqrt(3),(1-sqrt(3))/2,t,x');
        z = solveBurgers2(sqrt(3),(1+sqrt(3))/2,t,x');
        vRef = (w+z)/2;
        roRef = (z-vRef)/sqrt(3);
        pRef = roRef.^gamma;
        mRef = vRef.*roRef;
        ERef = pRef/(gamma-1) + 1/2*roRef.*vRef.^2;
    end

    %Error and convergence analysis
    if problem == 1 || problem == 2
        uRef = [roRef mRef ERef];

        if avoidShock == false
            uDif=uRef- u';
            L1Error(kk)= (1/(nSteps*3))*sum(abs(uDif),"all");
            L2Error(kk)= ((1/(nSteps*3))*sum(abs(uDif).^2,"all"))^(1/2);
            LinfError(kk) = max(max(abs(uDif)));
        end
        
        if avoidShock == true
            xShock = pi + alpha * tFinal;
            idx = find(abs(x-xShock)>1);
        
            uDif = uRef(idx) -  u(idx);
            L1Error(kk)= (1/numel(idx))*sum(abs(uDif));
            L2Error(kk)= (1/numel(idx))*sum(abs(uDif).^2)^(1/2);
            LinfError(kk) = max(abs(uDif));
        end
    end


    %Plotting at desired time
    if nSteps == nArray(1)
         fig_num = fig_num + 1;
         plot_components(u, uRef, x, kk, Method, tFinal,exact)
    end

    %Write to file

    % Open files for writing
fid1 = fopen('den.dat', 'w');
fid2 = fopen('pre.dat', 'w');
fid3 = fopen('vel.dat', 'w');

% Loop through data
for i = 1:nSteps
    % Write density
    fprintf(fid1, '%g %g\n', x(i), u(1,i));
    
    % Write pressure
    fprintf(fid2, '%g %g\n', x(i), (gamma - 1)*(u(3,i) - 0.5*u(2,i)^2/u(1,i)));
    
    % Write velocity
    fprintf(fid3, '%g %g\n', x(i), u(2,i)/u(1,i));
end

% Close files
fclose(fid1);
fclose(fid2);
fclose(fid3);

end

%Order Plot and Table
if problem ==1 || problem ==2
    fig_num = fig_num + 1;
    order_table_plot(Method, dx, tFinal, fig_num, nArray, order, L1Error, L2Error, LinfError)
end

function [fhatM, fhatP] = FD_WENO_Step(u,gam,nSteps,periodic)

        %Update variables with u
        [ro,m,E,v,p] = get_variables(u, gam);

        %Pad with periodicity
        if periodic == true
            ro = [ro(end - 2),ro(end - 1), ro(end), ro, ro(1),ro(2),ro(3)];
            m =  [m(end - 2),m(end - 1),  m(end),  m,  m(1), m(2), m(3)];
            E =  [E(end - 2),E(end - 1),  E(end),  E,  E(1), E(2), E(3)];
            v =  [v(end - 2),v(end - 1),  v(end),  v,  v(1), v(2), v(3)];
            p =  [p(end - 2),p(end - 1),  p(end),  p,  p(1), p(2), p(3)];
        elseif periodic == false %Change to set value from problem
            ro = [ro(1), ro(1), ro(1), ro, ro(end),ro(end),ro(end)];
            m =  [m(1),  m(1),  m(1),  m,  m(end), m(end), m(end)];
            E =  [E(1),  E(1),  E(1),  E,  E(end), E(end), E(end)];
            v =  [v(1),  v(1),  v(1),  v,  v(end), v(end), v(end)];
            p =  [p(1),  p(1),  p(1),  p,  p(end), p(end), p(end)];

        elseif periodic == 6 %Change to set value from problem
            ro = [1, 1, 1, ro, 1,1,1];
            m =  [0, 0,  0,  m,  0, 0, 0];
            E =  [1e-12,  1e-12,  1e-12,  E,  1e-12, 1e-12, 1e-12];
            %E =  [0,  0,  0,  E,  0, 0, 0];

            v =  [0,  0,  0,  v,  0, 0, 0];
            p =  [4e-13,  4e-13,  4e-13,  p,  4e-13, 4e-13, 4e-13];
            %p =  [0,  0,  0,  p,  0, 0, 0];

        end

        u = [ro;m;E];

        %Update variables with u
        fu = [m; (ro.*v.^2 + p); v.*(E + p)];

        %Build R and L using Roe Average 

        %Compute values for R
        t0 = sqrt(ro)./(sqrt(ro) + sqrt((circshift(ro,-1))));
        t1 = 1 - t0;
        Vxm = t0.*v + t1.*circshift(v,-1);

        H1 = (p + E)./ro;
        Hxm = t0.*H1 + t1.*circshift(H1,-1);

        qm = (0.5)*Vxm.^2;
        cm = sqrt((gam - 1)*(Hxm - qm));

        t2 = Vxm.*cm;
        
        %New alternative average
        % t0 = 0.5;
        % t1 = 0.5;
        % Vxm = v;
        % 
        % Hxm = (p + E)./ro;
        % %Hxm = t0.*H1 + t1.*circshift(H1,-1);
        % 
        % qm = (0.5)*Vxm.^2;
        % cm = sqrt((gam - 1)*(Hxm - qm));
        % 
        % t2 = Vxm.*cm;

        %End of alternative

        %Compute eigenvalues
        lam1 = Vxm - cm;
        lam2 = Vxm;
        lam3 = Vxm + cm;

        alpha1 = max(abs(lam1));
        alpha2 = max(abs(lam2));
        alpha3 = max(abs(lam3));

        %Local characteristic decomposition and LF
        fhatP = zeros(3,nSteps);
        for kk = 1:nSteps
            j = kk + 3;
            t2hold = t2(j); qmhold = qm(j);
            Vxmhold = Vxm(j); Hxmhold = Hxm(j); cmhold = cm(j); 
    
            R = zeros(3,3);
            R(1,:) = 1;
            R(2,1) = Vxmhold - cmhold;
            R(2,2) = Vxmhold;
            R(2,3) = Vxmhold + cmhold;
            R(3,1) = Hxmhold - t2hold;
            R(3,2) = qmhold;
            R(3,3) = Hxmhold + t2hold;
            
            rcmhold = 1/cmhold;
            b1hold = (gam - 1)*rcmhold^2;
            b2hold = qmhold*b1hold;
    
            t0hold = Vxmhold*rcmhold;
            t1hold = b1hold*Vxmhold;
            t2hold = 0.5*b1hold;
    
            L = zeros(3,3);
            L(1,1) = 0.5*(b2hold + t0hold);
            L(1,2) = -0.5*(t1hold + rcmhold);
            L(1,3) = t2hold;
            L(2,1) = 1 - b2hold;
            L(2,2) = t1hold;
            L(2,3) = -b1hold;
            L(3,1) = 0.5*(b2hold - t0hold);
            L(3,2) = -0.5*(t1hold-rcmhold);
            L(3,3) = t2hold;
            
    
            gp = zeros(5,3);
            gm = zeros(5,3);
            index = 0;
            for i = j-2:j+3
                index = index + 1;
    
                %LC
                w = L*u(:,i); %Account for padding in u
                g = L*fu(:,i);
    
                %LF
                %gp(index,:) = (1/2)*(g + [alpha1*w(1) ;alpha2*w(2); alpha3*w(3)]);
                gp(index,:) = (1/2)*(g + diag([alpha1, alpha2, alpha3])*w);
                gm(index,:) = (1/2)*(g - diag([alpha1, alpha2, alpha3])*w);
            end
    
            ghat = zeros(3,1);
    
            for l = 1:3
                [~,ghatp] = scalarWeno(gm(:,l));
                [ghatm,~] = scalarWeno(gp(:,l));
                ghat(l) = ghatp + ghatm;
            end
            
            fhatP(:,kk) = R*ghat;
        end

        fhatM = zeros(3,nSteps);
        
        if periodic == true
            for l = 1:3
                fhatM(l,:) = circshift(fhatP(l,:),1);
            end       
        else 
             %Local characteristic decomposition and LF
             for l = 1:3
                fhatM(l,2:end) = fhatP(l,1:end-1);
             end  

        for kk = 1
            j = kk + 2;
            t2hold = t2(j); qmhold = qm(j);
            Vxmhold = Vxm(j); Hxmhold = Hxm(j); cmhold = cm(j); 
    
            R = zeros(3,3);
            R(1,:) = 1;
            R(2,1) = Vxmhold - cmhold;
            R(2,2) = Vxmhold;
            R(2,3) = Vxmhold + cmhold;
            R(3,1) = Hxmhold - t2hold;
            R(3,2) = qmhold;
            R(3,3) = Hxmhold + t2hold;
            
            rcmhold = 1/cmhold;
            b1hold = (gam - 1)*rcmhold^2;
            b2hold = qmhold*b1hold;
    
            t0hold = Vxmhold*rcmhold;
            t1hold = b1hold*Vxmhold;
            t2hold = 0.5*b1hold;
    
            L = zeros(3,3);
            L(1,1) = 0.5*(b2hold + t0hold);
            L(1,2) = -0.5*(t1hold + rcmhold);
            L(1,3) = t2hold;
            L(2,1) = 1 - b2hold;
            L(2,2) = t1hold;
            L(2,3) = -b1hold;
            L(3,1) = 0.5*(b2hold - t0hold);
            L(3,2) = -0.5*(t1hold-rcmhold);
            L(3,3) = t2hold;
            
    
            gp = zeros(5,3);
            gm = zeros(5,3);
            index = 0;
            for i = j-2:j+3
                index = index + 1;
    
                %LC
                w = L*u(:,i); %Account for padding in u
                g = L*fu(:,i);
    
                %LF
                %gp(index,:) = (1/2)*(g + [alpha1*w(1) ;alpha2*w(2); alpha3*w(3)]);
                gp(index,:) = (1/2)*(g + diag([alpha1, alpha2, alpha3])*w);
                gm(index,:) = (1/2)*(g - diag([alpha1, alpha2, alpha3])*w);
            end
    
            ghat = zeros(3,1);
    
            for l = 1:3
                [~,ghatp] = scalarWeno(gm(:,l));
                [ghatm,~] = scalarWeno(gp(:,l));
                ghat(l) = ghatp + ghatm;
            end
            
            fhatM(:,1) = R*ghat;
        end
  

        end

end

function [umOut,upOut] = scalarWeno(u)
    umm = u(1);
    um = u(2);
    up = u(4);
    upp = u(5);
    uppp = u(6);
    u = u(3);

    um1 = 1/3*umm - 7/6*um + 11/6*u; 
    um2 = -1/6*um + 5/6*u + 1/3*up; 
    um3 = 1/3*u + 5/6*up - 1/6*upp; 
    
    %Get proper coefficients
    c(1) = 1/10;
    c(2) = 3/5;
    c(3) = 3/10;
    
   %Get proper coefficients
    beta1 = 13/12*( umm - 2*um + u ).^2 + 1/4*( umm - 4*um + 3*u ).^2; 
    beta2 = 13/12*( um - 2*u + up ).^2 + 1/4*( um - up ).^2; 
    beta3 = 13/12*( u - 2*up + upp ).^2 + 1/4*( 3*u - 4*up + upp ).^2;
    
    %Get w tilde and normalize
    eps = 10e-6;
    wt1 = c(1) ./ (beta1 + eps).^2;
    wt2 = c(2) ./ (beta2 + eps).^2;
    wt3 = c(3) ./ (beta3 + eps).^2;
    wt = [wt1, wt2, wt3];
    w1 = wt(:,1) ./ sum(wt, 2);
    w2 = wt(:,2) ./ sum(wt, 2);
    w3 = wt(:,3) ./ sum(wt, 2);
    
    %Final weight stencil for um 
    
    umOut = w1.*um1 + w2.*um2 + w3.*um3; 
    
    %Use symmetry

    up1 = 11/6*up - 7/6*upp + 1/3*uppp; 
    up2 = 1/3*u + 5/6*up - 1/6*upp; 
    up3 = -1/6*um + 5/6*u + 1/3*up;
    
    %Smoothness indicators
    beta1 = 13/12*( uppp - 2*upp + up ).^2 + 1/4*( uppp - 4*upp + 3*up ).^2; 
    beta2 = 13/12*( upp - 2*up + u ).^2 + 1/4*( upp - u ).^2; 
    beta3 = 13/12*( up - 2*u + um ).^2 + 1/4*( 3*up - 4*u + um ).^2;
    
    %Get w tilde and normalize
    eps = 10e-6;
    wt1 = c(1)./(beta1 + eps).^2;
    wt2 = c(2)./(beta2 + eps).^2;
    wt3 = c(3)./(beta3 + eps).^2;
    wt = [wt1, wt2, wt3];
    w1 = wt(:,1) ./ sum(wt, 2);
    w2 = wt(:,2) ./ sum(wt, 2);
    w3 = wt(:,3) ./ sum(wt, 2);
    
    %Final weight stencil for up
    upOut = w1.*up1 + w2.*up2 + w3.*up3; 
end



