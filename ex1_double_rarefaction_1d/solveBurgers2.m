function u = solveBurgers2(alpha, beta, t, x)

%Using change of variables, method of
%characteristics and Newton's iteration to solve Burger's equation with
%general periodic initial condition

 %Locate the shock
    if (pi+alpha*t>2*pi)
        x = x + 2*pi*floor((pi+alpha*t)/(2*pi));
    end

 %%%%Change of variables and move into (-pi,pi) interval%%%%%

    
    dist = pi+alpha*t;
    %Bring the x_1 into the interval (-pi,pi)
    x_1= x(x<dist);
    
    if (~isempty(x(x<dist)))
        %Transform into v coordinates (x_1, t_1) -> (x'_1, t'_1)
        x1_1 = x_1-alpha*t;
    
        %Solve Burger's in B coordinates
        u = burgers(beta,x1_1,t);
        v = alpha+u;
            
        if (~isempty(x(x>dist)))
            x_2 = x(x>=dist);
            x_2 = flip(2*dist-x_2);
            x1_2 = x_2-alpha*t;
    
            u2 = burgers(beta,x1_2,t);
            v2 = alpha+u2;
    
            u = [v; flip(2*alpha-v2)];
        else
            u = v;
        end
    elseif (~isempty(x(x>dist)))
        x_2 = x(x>=dist);
        x_2 = flip(2*dist-x_2);
        x1_2 = x_2-alpha*t;
    
        u2 = burgers(beta,x1_2,t);
        v2 = alpha+u2;

        u = flip(2*alpha-v2);
    else
        v=0;
        u = alpha + beta*v;
        return        
    end
    
    
    return
end

function [f,fp] = myfun(xx,x11,t,b)

f = x11 - xx - (b*sin(xx)).*t;
fp = -1 - t.*b.*cos(xx);

return
end

function u = burgers(b,x,t)

% compute x0 using nonlinear solve
x_next = zeros(length(x),1);

for j=1:length(x)
    x_next(j) = fzero(@(xx)myfun(xx,x(j),t,b),0);
end

u = b*sin(x_next);

return

end