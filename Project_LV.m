%% point 2.1
clc
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MARKET DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Spot price
S_0=3583.67;

%Market expires
T=[ 0.090 	 0.247 	 0.496 	 0.764 	 1.014 	 1.510 	 2.008 	 2.507];

% forwards at market expiries
Fwd = [3581.24	3573.30	3490.39	3480.82	3463.77	3369.32	3346.41	3255.06];

% discount factor at market expiry
DF = [ 1.0003610 	 1.0008430 	 1.0013930 	 1.0021070 	 1.0028740 	 1.0041040 	 1.0046760 	 1.0039930];

% market strikes
K = [
    3300	3000	2600	2350	2200	1900	1750	1550;
    3500	3400	3200	3100	3050	2850	2750	2600;
    3550	3500	3400	3350	3300	3200	3150	3000;
    3600	3550	3500	3500	3450	3350	3350	3250;
    3650	3650	3550	3600	3600	3550	3550	3500;
    3700	3700	3700	3750	3800	3800	3850	3850;
    3750	3900	3950	4100	4200	4400	4600	4800];

% market volatilities
MktVol = [0.1652	0.2181	0.2549	0.2744	0.2749	0.2839	0.2814	0.288;
          0.1074	0.1509	0.177	0.1889	0.1924	0.2028	0.2043	0.2073;
          0.0935	0.1351	0.1529	0.1634	0.1716	0.1788	0.1822	0.1877;
          0.0832	0.1274	0.1412	0.1488	0.1597	0.1691	0.1723	0.1773;
          0.0777	0.1137	0.1357	0.1401	0.1486	0.1573	0.1634	0.1682;
          0.0768	0.1078	0.1223	0.1295	0.1364	0.1458	0.1527	0.1587;
          0.0793	0.0984	0.1085	0.1144	0.1199	0.1309	0.1375	0.1475];

   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALIBRATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calibrate model parameters
[r,q] = calibrate_r_q(S_0,T,DF,Fwd);

% normalized market strikes
[rows, cols] = size(K);
K_norm = K ./ repmat(Fwd, rows, 1);

% Dupire solver settings
Lt = 10;
Lh = 500;
K_min = 0.1;
K_max = 3;
Scheme = 'cn';

% calibration settings
Threshold = 0.0010;
MaxIter = 200;

[V, ModelVol, MaxErr] = calibrator(T,K_norm,MktVol,Threshold,MaxIter,Lt,Lh,K_min,K_max,Scheme);

% plot local volatility function vs market implied volatility
figure;
plot(K(:,1),MktVol(:,1),'o',K(:,1),ModelVol(:,1),':.',K(:,1),V(:,1),':.b','linewidth',1.5);
xline(S_0,'-',{'Spot price'});
title('Calibrated model and local volatility for asset E CORP');
legend('MktVol','ModelVol','LocalVol')

figure;
plot(MaxErr,'.','MarkerSize',15);
title('calibration error at each iteration of the fixed-point calibration');

%% point 2.2
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRICING CALL OPTIONS THROUGH MC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

S_0=3583.67;
k=[0.9 1.1];

lv_price=[];
lv_impl_vol=[];
price_black=[];
black_impl_vol=[];

for i=1:length(k)
    
expiry=0.5;
strike=k(i)*S_0;

N = 1000000; %MC simulations
M = 50; %timesteps

% MC simulation (Local Vol)
S = lv_simulation_log(T,S_0,r,q,V,K,N,M,expiry);

% LV price of a call option
discount_factor = interp1(T,DF,expiry);
lv_price = [lv_price,discount_factor*mean(max(S(1,:) - strike,0))];
fprintf('Price of the call with perc_strike equal to %.2f\n', k(i));
disp(lv_price(i));



% LV implied volatility
fwd = interp1(T,Fwd,expiry);
lv_impl_vol = [lv_impl_vol, blsimpv(fwd,strike,0,expiry,lv_price/discount_factor)];
fprintf('Implied volatility of the call with perc_strike equal to %.2f\n', k(i));
disp(lv_impl_vol(i));

%{
% MC simulation (Black)
time_idx = find(T>=expiry,1);
sigma = interp1(K(:,time_idx),MktVol(:,time_idx),strike); % model parameter
N_B = 1000000; %MC simulations
S_B = black_simulation_log(T,Fwd,sigma,N_B,expiry); 

% Black price of a call option
price_black = [price_black,discount_factor*mean(max(S_B(1,:) - strike,0))];
black_impl_vol = [black_impl_vol, blsimpv(fwd,strike,0,expiry,price_black/discount_factor)];


% Probability densities generated by Monte Carlo simulation for LV and
% Black model
figure;
hist(S_B(1,:),100);
title('Probability density generated by Black model');

figure;
hist(S(1,:),100);
title('Probability density generated by LV model') 
%}
end
%% point 2.3
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRICING FOWARDS STARTING OPTIONS THROUGH MC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
P=[];
model_impl_vol=[];

for i=1:length(k)
% option data
expiry = [2 2.5];

% MC settings
N = 1000000; %MC simulations
M = 100; %timesteps

% MC simulation (LV)
S = lv_simulation_log(T,S_0,r,q,V,K,N,M,expiry);

% option price
discount_factor = discount(T,DF,expiry(2));

P = [P,discount_factor*mean(max(S(2,:) - k(i)*S(1,:),0))];
fprintf('Price of the fs with perc_strike equal to %.2f\n', k(i));
disp(P(i));

fwd(1) = forward(S_0,T,r,q,expiry(1));
fwd(2) = forward(S_0,T,r,q,expiry(2));
model_impl_vol = [model_impl_vol, blsimpv(fwd(2),k(i)*fwd(1),0,expiry(2)-expiry(1),P/discount_factor)];
fprintf('Implied volatility of the fs with perc_strike equal to %.2f\n', k(i));
disp(model_impl_vol(i));

end

%% point 2.4
clc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRICING SPOT START CALL OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
% spot start option data
Expiry = 0.5;

% compute price of spot-start call options 
[ k, C ] = solve_dupire( T, K_norm, V, Expiry, Lt, Lh, K_min, K_max, Scheme);    

% compute model implied volatilities
perc_strikes_spot_start=[];
model_impl_vol_spot_start=[];
for i=1:length(k)
    if k(i)>0.9001 && k(i)<1.1
        perc_strikes_spot_start = [perc_strikes_spot_start k(i)];
        model_impl_vol_spot_start = [model_impl_vol_spot_start blsimpv(1,k(i),0,Expiry,C(i))];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRICING FORWARD START CALL OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% option prices
perc_strikes = 0.9:0.01:1.1;
model_impl_vol=[];
for x = perc_strikes
    P = discount_factor*mean(max(S(2,:) - x*S(1,:),0));
    model_impl_vol = [model_impl_vol, blsimpv(fwd(2),x*fwd(1),0,expiry(2)-expiry(1),P/discount_factor)];
end

plot(perc_strikes_spot_start,model_impl_vol_spot_start,perc_strikes,model_impl_vol, 'linewidth',1.5);
title('Model implied vol with option maturity 6m');
legend('Spot impl vol','Fwd impl vol');

fprintf('skewness spot impl vol')
skewness(model_impl_vol_spot_start)
fprintf('skewness fwd impl vol')
skewness(model_impl_vol)

%% point 3
clc 
clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MARKET DATA (FAIL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%market expiries
T_fail = [ 0.01644 	 0.03288 ];

% forwards at market expiries
Fwd_fail = [4.50708	4.50708];

% discount factor at market expiries
DF_fail = [ 1.0000000 	 1.0000000 ];

% market strikes
K_fail = [ 4.41060	4.41060;
           4.45479	4.45479;
           4.50809	4.50809;
           4.56930	4.56930;
           4.63267	4.63267 ];

% market volatilities
MktVol_fail = [0.1903	0.1203;
               0.1975	0.1275;
               0.2050	0.1350;
               0.4175	0.1475;
               0.2323	0.1623];
% spot
S_0 = 4.50708;

% calibrate model parameters
[r,q] = calibrate_r_q(S_0,T_fail,DF_fail,Fwd_fail);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALIBRATION (FAIL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% normalized market strikes
[rows, cols] = size(K_fail);
K_norm = K_fail ./ repmat(Fwd_fail, rows, 1);

% Dupire solver settings
Lt = 10;
Lh = 200;
K_min = 0.1;
K_max = 3;
Scheme = 'cn';

% calibration settings
Threshold = 0.0010;
MaxIter = 100;

% [V, ModelVol, MaxErr] = calibrator(T_fail,K_norm,MktVol_fail,Threshold,MaxIter,Lt,Lh,K_min,K_max,Scheme);

price_t1 = [];
price_t2 = [];
for i = 1:5

price_t1(i) = blsprice(Fwd_fail(1),K_fail(i,1),0,T_fail(1),MktVol_fail(i,1));
price_t2(i) = blsprice(Fwd_fail(2),K_fail(i,2),0,T_fail(2),MktVol_fail(i,2));

end

figure;
hold on
plot(K_fail(:,1),price_t1,'o-','LineWidth',2)
plot(K_fail(:,1),price_t2,'o-','LineWidth',2)
title('convexity of the map Ki,j → C_0(Ti,Ki,j) for fixed T_i')
legend('Prices for T=T_1','Prices for T=T_2')


%% point 4.1

clc
clearvars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MARKET DATA (FX)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Market expiries
T = [ 0.047 	 0.093 	 0.178 	 0.258 	 0.510 	 0.762 	 1.011 ];

%domestic discount factor
disc_fact=[ 0.9993910 	 0.9986680 	 0.9972880 	 0.9958860 	 0.9910190 	 0.9855220 	 0.9804450];

%foreign discount factor
disc_fact_for=[ 1.0002703 	 1.0003940 	 1.0006179 	 1.0008642 	 1.0015410 	 1.0018667 	 1.0026841];

% forwards at market expiries
Fwd = [1.181939	 1.182941	1.184843	1.186803	1.193438	1.200485	1.207686];

% spot rate
spot = 1.18090;

% market deltas
Delta = [ 0.1 0.25 0.5 0.75 0.9];

% market volatilities
MktVol = [0.072842	 0.064896	0.065436	0.076120	0.075762	0.077481	0.079380;
           0.067591	 0.057149	0.061758	0.072632	0.072016	0.074099	0.075251;
           0.064659	 0.057379	0.060844	0.070767	0.070899	0.073041	0.074308;
           0.061035	 0.059441	0.063119	0.072028	0.073410	0.075848	0.077449;
           0.067613	 0.061692	0.066479	0.075432	0.078536	0.080742	0.083532];

% market strikes
K = zeros(length(Delta),length(T));

% find K such that BS_Delta(K,Fwt,T,MktVol) = Delta
for i = 1:length(Delta)
    for j = 1:length(T)
        K(i,j) = fzero(@(Strike) blsdelta(Fwd(j),Strike,0,T(j),MktVol(i,j))-(1-Delta(i)), Fwd(j));
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALIBRATION (FX)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% calibrate model parameters
[r,q] = calibrate_r_q(spot,T,disc_fact,Fwd);

% normalized market strikes
[rows, cols] = size(K);
K_norm = K ./ repmat(Fwd, rows, 1);

% Dupire solver settings
Lt = 10;
Lh = 500;
K_min = 0.5;
K_max = 2.5;
Scheme = 'cn';

% calibration settings
Threshold = 0.001;
MaxIter = 100;

[V, ModelVol, MaxErr] = calibrator(T,K_norm,MktVol,Threshold,MaxIter,Lt,Lh,K_min,K_max,Scheme);

% plot local volatility function vs market implied volatility
figure;
plot(K(:,1),MktVol(:,1),'o',K(:,1),ModelVol(:,1),':.',K(:,1),V(:,1),':.b','linewidth',1.5);
xline(spot,'-',{'Spot Rate'});
title('Calibrated model and local volatility for asset EUR/USD');
legend('MktVol','ModelVol','LocalVol');

figure;
plot(MaxErr,'.','MarkerSize',15);
title('calibration error at each iteration of the fixed-point calibration');

%% point 4.2
clc

%Let's consider a plain vanilla option and a digital option with the following caratheristics:
expiry=T(5);
strike=K(5,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRICING (PLAIN VANILLA OPTION)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MC settings
N = 100000; %MC simulations
M = 100; %timesteps

% MC simulation (LV)
S = lv_simulation_log(T,spot,r,q,V,K,N,M,expiry);

% option price (LV)
discount_factor = interp1(T,disc_fact,expiry);
lv_price_pv = discount_factor*mean( max(S(1,:) - strike, 0) );

fprintf('Price of the plain vanilla option with strike=%.2f and expiry=%.1f\n',strike,expiry);
disp(lv_price_pv);

%Confidence interval of 95 percent
z=1.96;  %quantile 0.95
v_n_pv=std(discount_factor*max(S(1,:) - strike, 0)); %sample variance of the Monte Carlo simulation
Error_pv=z*sqrt(v_n_pv)/sqrt(N);

fprintf('Confidence interval of level 95 for the option price\n');
fprintf('%.5f $ +/- %.5f $\n\n',lv_price_pv,Error_pv);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRICING (DIGITAL OPTION)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

lv_price_do= discount_factor*mean((S(1,:) > strike));

fprintf('Price of the digital option with strike=%.2f and expiry=%.1f\n',strike,expiry);
disp(lv_price_do);

%Confidence interval of 95 percent
z=1.96;  %quantile 0.95
v_n_do=std(discount_factor*(S(1,:) > strike)); %sample variance of the Monte Carlo simulation
Error_do=z*sqrt(v_n_do)/sqrt(N);

fprintf('Confidence interval of level 95 for the digital option price\n');
fprintf('%.5f $ +/- %.5f $\n\n',lv_price_do,Error_do);

%% point 4.3
clc

sigma = MktVol(5,2);
expiry=T(5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRICING UNDER BLACK (PLAIN VANILLA OPTION)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MC simulation (Black)
N_B = 100000; %MC simulations
S_B = black_simulation_log(T,Fwd,sigma,N_B,expiry); 

% Black price of a plain vanilla option
price_black_pv = discount_factor*mean(max(S_B(1,:) - strike,0));
%black_impl_vol = blsimpv(Fwd,strike,0,expiry,price_black/discount_factor);

fprintf('Price (under BS) of the plain vanilla option with strike=%.2f and expiry=%.1f\n',strike,expiry);
disp(price_black_pv);

%Confidence interval of 95 percent
z=1.96;  %quantile 0.95
v_black_pv=std(discount_factor*max(S_B(1,:) - strike, 0)); %sample variance of the Monte Carlo simulation
Error_black_pv=z*sqrt(v_black_pv)/sqrt(N_B);

fprintf('Confidence interval of level 95 for the option price (under BS):\n');
fprintf('%.5f $ +/- %.5f $\n\n',price_black_pv,Error_black_pv);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRICING UNDER BLACK(DIGITAL OPTION)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

price_black_do= discount_factor*mean((S_B(1,:) > strike));

fprintf('Price (under BS) of the digital option with strike=%.2f and expiry=%.1f\n',strike,expiry);
disp(price_black_do);

%Confidence interval of 95 percent
z=1.96;  %quantile 0.95
v_black_do=std(discount_factor*(S_B(1,:) > strike)); %sample variance of the Monte Carlo simulation
Error_black_do=z*sqrt(v_black_do)/sqrt(N_B);

fprintf('Confidence interval of level 95 for the digital option price\n');
fprintf('%.5f $ +/- %.5f $\n\n',price_black_do,Error_black_do);