function X = simulate_closed_loop(Fcl, x0, T)
    nx = size(Fcl,1);
    X = zeros(nx, T);
    X(:,1) = x0;
    for k = 2:T
        X(:,k) = Fcl * X(:,k-1);
    end
end
