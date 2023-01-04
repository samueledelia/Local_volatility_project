function [V, ModelVol, MaxErr] = calibrate_fp(T,K_norm,MktVol,threshold,MaxIter,Lt,Lh,K_min,K_max,Scheme)

% calibrates a local volatility model using fixed-point algorithm. 
%Code is optimized so that at each iteration of the algorithm only 
%one Dupire equation is solved

nb_iter = 1;
MaxErr(nb_iter) = 100;

% initial guess for the LV matrix
V = MktVol;
% initialize ModelVol just once at the beginning
[rows, cols] = size(K_norm);
ModelVol = zeros(rows,cols);

% forward market volatility
mkt_fwd_vol = fwd_from_spot_vol(T,MktVol);

% calibration progress
h = waitbar(0,'Calibration...');

idx = 1;
MaxIter = MaxIter*length(T);


while ( MaxErr(nb_iter) > threshold && nb_iter < MaxIter )

    % update iteration number
    % MaxIter is referred to the orginal meaning of it. Actual iterations
    % are length(T) times more
    nb_iter = nb_iter+1;
    
    % compute model implied volatilities with one Dupire
    expiry = T(idx); %solve Dupire equation only with expiry of one Ti
    [ k, C ] = solve_dupire( T, K_norm, V, expiry, Lt*idx, Lh, K_min, K_max, Scheme);
    price = interp1(k,C,K_norm(:,idx),'linear');
    ModelVol(:,idx) = blsimpv(1,K_norm(:,idx),0,expiry,price); %only one column is modified
    
    % compute max_err
    MaxErr(nb_iter) = max(max(abs(ModelVol-MktVol))); 
    %(^the error in the first length(T) iteration will be huge)
    waitbar(threshold/MaxErr(nb_iter));
    
    % compute model fwd vol 
    model_fwd_vol = fwd_from_spot_vol(T,ModelVol);
    
    % compute new LV parameters
    V = V ./ model_fwd_vol .* mkt_fwd_vol;
    
    disp('fine')
    disp(nb_iter-1)
    
    if idx == length(T)
        idx = 1;
    else
        idx = idx+1;
    end

end

close(h);
MaxErr = MaxErr(2:nb_iter);

end