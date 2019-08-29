function [l, u] = index_lower_upper_bounds(om)

    %% number of variables
    nlin = om.lin.N;        %% number of linear constraints
    
    u = Inf(nlin, 1);       %% upper bound
    l = -u;                 %% lower bound
    
    %% fill in each piece
    for k = 1:om.lin.NS
        name = om.lin.order(k).name;
        idx  = om.lin.order(k).idx;
        [~, lk, uk, ~, i1, iN] = retrieve_lin_constraints(om, name, idx);
        
%         l(i1:iN) = lk;
%         u(i1:iN) = uk;

          l(i1:iN) = k;
          u(i1:iN) = k;
        
    end
    
end