function[S] = lv_simulation_log(T,Spot,r,q,V,K,N,M,expiry)
   % Monte Carlo simulation for the asset S under the LV model
   % T.. LV expiries
   % Spot.. spot of the asset needed to reconstruct forwards
   % r.. risk-free rates at T, q.. dividend yiels at T
   % V.. LV matrix, K.. LV strikes
   % M.. MC timesteps, N.. MC simulations, expiry.. MC time horizon

   logX = zeros(N,1);
   S = zeros(length(expiry),N);

   % model forwards at T
   Fwd = zeros(1,length(T));
   for i=1:length(T)
      Fwd(i) = forward(Spot,T,r,q,T(i));
   end
   
   % normalized LV strikes
   [rows, ~] = size(K);
   K_norm = K ./ repmat(Fwd, rows, 1);

   % generate random numbers
   W = randn(N,M,length(expiry));

   for k = 1:length(expiry)
      if k==1
          dt = expiry(k)/M;
          t = (0:M)'*dt;
      else
          dt = (expiry(k)-expiry(k-1))/M;
          t = expiry(k-1) + (1:M)'*dt;
      end

      for j = 1:M
          eta = localvol(T,K_norm,V,t(j),exp(logX(:)));
          logX(:) = logX(:) - 0.5 * dt * eta.^2 + eta .* W(:,j,k) * sqrt(dt);
      end

      fwd = forward(Spot,T,r,q,expiry(k));
      S(k,:) = fwd * exp(logX(:));
   end
end

