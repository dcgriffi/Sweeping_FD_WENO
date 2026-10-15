function order_table_plot(Method, dx, tFinal, fig_num, nArray, order, L1Error, L2Error, LinfError)
    p = order;
    figure(fig_num)
    hold on;
    plot(log10(dx),log10(L2Error),'.-','LineWidth',2,'MarkerSize',20)
    plot(log10(dx),p*log10(dx),'k','LineWidth',3) %Reference Slope = 3
    xlabel('Log_{10}(\Deltax)')
    ylabel('Log_{10}(Error)')
    title("L2 Convergence Plot for " + num2str(Method) + " at time " +tFinal)
    legend('', 'p=3','location','southeast')
    
    L1Order = (log10(L1Error(1:end-1))-log10(L1Error(2:end)))./(log10(dx(1:end-1))-log10(dx(2:end))) ;
    L2Order = (log10(L2Error(1:end-1))-log10(L2Error(2:end)))./(log10(dx(1:end-1))-log10(dx(2:end))) ;
    LinfOrder = (log10(LinfError(1:end-1))-log10(LinfError(2:end)))./(log10(dx(1:end-1))-log10(dx(2:end))); 
    
    %Format table
    fractions = arrayfun(@(x) sprintf('1/%d', x), nArray, 'UniformOutput', false);
    T = table(fractions', L1Error(:),[0;L1Order(:)], L2Error(:),[0;L2Order(:)], LinfError(:),[0;LinfOrder(:)], ...
        'VariableNames',{'dx','L1 Error','L1 Order','L2 Error','L2 Order','Linf Error','Linf Order'});
    disp(T)
end