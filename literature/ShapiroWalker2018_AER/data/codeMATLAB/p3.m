pm.alpha_nj    = repmat(pm.alpha',[N 1]);
pm.sigma_nj    = repmat(pm.sigma',[N 1]);
pm.theta_nj    = repmat(pm.theta',[N 1]);
pm.alpha_jy    = repmat(pm.alpha ,[1 Y]);
pm.sigma_jy    = repmat(pm.sigma ,[1 Y]);
pm.theta_jy    = repmat(pm.theta ,[1 Y]);
pm.alpha_njy   = repmat(pm.alpha',[N 1 Y]);
pm.sigma_njy   = repmat(pm.sigma',[N 1 Y]);
pm.theta_njy   = repmat(pm.theta',[N 1 Y]);
pm.alpha_nnj   = repmat(reshape(pm.alpha,[1 1 J]),[N N 1]);
pm.sigma_nnj   = repmat(reshape(pm.sigma,[1 1 J]),[N N 1]);
pm.theta_nnj   = repmat(reshape(pm.theta,[1 1 J]),[N N 1]);
pm.alpha_nnjy  = repmat(reshape(pm.alpha,[1 1 J 1]),[N N 1 Y]);
pm.sigma_nnjy  = repmat(reshape(pm.sigma,[1 1 J 1]),[N N 1 Y]);
pm.theta_nnjy  = repmat(reshape(pm.theta,[1 1 J 1]),[N N 1 Y]);


Eds         = squeeze(sum(xbilat,1));           
Rds         = squeeze(sum(xbilat,2));           
beta        = Eds./repmat(sum(Eds,2),[1 J 1]);  


lambda = xbilat ./ repmat(sum(xbilat,1),[N 1 1 1]);
zeta   = xbilat ./ repmat(sum(xbilat,2),[1 N 1 1]); 
Ed     = squeeze(sum(sum(xbilat,3),1));
Rd     = squeeze(sum(sum(xbilat,3),2));
NXds   = Rds - Eds;                                 
NXd    = Rd - Ed;                                   
NXAds  = squeeze(sum(NXds.*(pm.sigma_njy-1).*(pm.theta_njy-pm.alpha_njy+1) ./(pm.theta_njy.*pm.sigma_njy),2));    


for j = 1:J
    Z_co2(j,:) = interp1([1991,1994,1998,2002,2006,2010]',co2_poll(j,:)',(1990:2008)','linear','extrap')';
for n = 1:length(pollutants)
    Z(j,:,n) = interp1q([1990,1996,1999,2002,2005,2008]',poll(j,:,n)',(1990:2008)')';
end
end
Z(:,:,7) = Z_co2;

    
bl.xbilat0  = xbilat(:,:,:,1);
bl.lambda0  = lambda(:,:,:,1);
bl.zeta0    = zeta(:,:,:,1);
bl.beta0    = beta(:,:,1);
bl.Rd0      = Rd(:,1);
bl.Rds0     = Rds(:,:,1);
bl.NXd0     = NXd(:,1);
bl.NXAds0   = NXAds(:,1);
bl.wwM_hat0 = ones(N*J+N-1,1);

w_hat           = squeeze(sum(vship,2)) ./ repmat(sum(vship(:,:,1),2),[1 Y]);
lambda_hat      = lambda                ./ repmat(lambda(:,:,:,1),[1 1 1 Y]);
shocks.beta_hat = beta                  ./ repmat(beta(:,:,1)    ,[1 1   Y]);
Z_hat           = Z                     ./ repmat(Z(:,1,:),[1 Y 1]);
Rds_hat         = Rds                   ./ repmat(Rds(:,:,1),[1 1 Y]);
shocks.NXd_hat  = NXd                   ./ repmat(NXd(:,1),[1 Y]); 
shocks.NXAds_hat = NXAds                ./ repmat(NXAds(:,1),[1 Y]);


w_oNJY= repmat(permute(w_hat,[1 3 2]),[1 J 1]);
M_hat = Rds_hat ./ w_oNJY;


pwrE          = 1-pm.theta_nnjy./((pm.sigma_nnjy-1).*(1-pm.alpha_nnjy));
M_hat_nnjy    = repmat(permute(M_hat   ,[1 4 2 3]),[1 N 1 1]);
w_hat_nnjy    = repmat(permute(w_hat   ,[1 3 4 2]),[1 N J 1]);
w_hatd_nnjy   = repmat(permute(w_hat   ,[3 1 4 2]),[N 1 J 1]);
beta_hat_nnjy = repmat(permute(shocks.beta_hat,[4 1 2 3]),[N 1 1 1]);
Rd_nnjy       = repmat(permute(Rd      ,[3 1 4 2]),[N 1 J 1]);
NXd_nnjy      = repmat(permute(NXd     ,[3 1 4 2]),[N 1 J 1]);
bl.Rd0_nnjy   = repmat(        bl.Rd0'            ,[N 1 J Y]);
bl.NXd0_nnjy  = repmat(        bl.NXd0'           ,[N 1 J Y]);
shocks.Gamma_hat_star = (lambda_hat./ (M_hat_nnjy .* (w_hat_nnjy.^(-pm.theta_nnjy)))) .* (((beta_hat_nnjy./w_hatd_nnjy) .*(Rd_nnjy-NXd_nnjy)./(bl.Rd0_nnjy-bl.NXd0_nnjy)) .^ pwrE);
t_hat = squeeze(M_hat(us,:,:)) .* repmat(w_hat(us,:),[J 1]) ./ Z_hat(:,:,mainpoll);
shocks.Gamma_hat_t = ones(size(shocks.Gamma_hat_star));
shocks.Gamma_hat_t(us,:,:,:) = repmat(permute( t_hat.^(-pm.alpha_jy.*pm.theta_jy./(1-pm.alpha_jy)) ,[3 4 1 2]),[1 N 1 1]);

shocks.Gamma_hat_foreign  = [shocks.Gamma_hat_star(1:(us-1),:,:,:);ones(1,N,J,Y)];
shocks.Gamma_hat_domestic = [ones(1:(us-1),N,J,Y);shocks.Gamma_hat_star(us,:,:,:)./shocks.Gamma_hat_t(us,:,:,:)];
    
for n = 1:Y
    wwM_hat(:,n) = [w_hat(2:end,n);reshape(M_hat(:,:,n),[N*J 1])];
end

if     mainpoll == 1 pollname = 'co';
elseif mainpoll == 2 pollname = 'nox';
elseif mainpoll == 3 pollname = 'pm10';
elseif mainpoll == 4 pollname = 'pm25';
elseif mainpoll == 5 pollname = 'so2';
elseif mainpoll == 6 pollname = 'voc';
elseif mainpoll == 7 pollname = 'co2';
end

save temp.mat


for n = 1:Y
for loop_shock = 1:5
    
    [mywwM_hat,fval,exitflag] = fsolve(@(wwM_hat)p4(wwM_hat, bl, shocks, loop_shock, n, N, J, pm),wwM_hat(:,1), options);
    resvecF = sum(fval.^2);

    ww      = mywwM_hat(1:(N-1));
    M_hat   = reshape(mywwM_hat(N:end),[N J]);
    wA      = (sum(bl.Rd0)./bl.Rd0(1)) - sum(ww.*bl.Rd0(2:end)./bl.Rd0(1)); 
    w_hat   = [wA;ww];

    if loop_shock == 3; Z_hat_cf(:,n,loop_shock) = (w_hat(us).*M_hat(us,:)'./t_hat(:,n)); end;
    if loop_shock ~= 3; Z_hat_cf(:,n,loop_shock) = (w_hat(us).*M_hat(us,:)'./t_hat(:,1)); end;

    clear Gamma_hat_s
end
end
clearvars -except Z_hat_o_cf Z_hat_cf minimized_fval resvecF
load temp.mat           

bl.Z0 = repmat(Z(:,1,:),[1 Y 1]);
Z_hat_o_cf = 100 .* squeeze( sum(Z_hat_cf .* repmat(bl.Z0(:,:,mainpoll),[1 1 5]) ,1) ./ sum(repmat(bl.Z0(:,:,mainpoll),[1 1 5]),1)); 
Z_hat_o = squeeze(100.*sum(Z,1)./sum(bl.Z0,1));



