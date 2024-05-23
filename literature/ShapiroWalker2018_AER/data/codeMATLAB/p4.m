function q = p4( wwM_hat, bl, shocks, loop_shock, n, N, J, pm)           

    NXd_hat    = ones(N,1);
    NXAds_hat   = ones(N,1);
    Gamma_hat  = ones(N,N,J);
    beta_hat   = ones(N,J);
    
    if loop_shock == 1 Gamma_hat  = shocks.Gamma_hat_foreign(:,:,:,n); end;
    if loop_shock == 2 Gamma_hat  = shocks.Gamma_hat_domestic(:,:,:,n); end;
    if loop_shock == 3 Gamma_hat  = shocks.Gamma_hat_t(:,:,:,n);       end;
    if loop_shock == 4 beta_hat   = shocks.beta_hat(:,:,n);            end;
    if loop_shock == 5 NXd_hat    = shocks.NXd_hat(:,n);               end;
    if loop_shock == 5 NXAds_hat   = shocks.NXAds_hat(:,n);            end;
    
    if loop_shock == 6 Gamma_hat  = shocks.Gamma_hat_star(:,:,:,n);    end;
    if loop_shock == 6 beta_hat   = shocks.beta_hat(:,:,n);            end;
    if loop_shock == 6 NXd_hat    = shocks.NXd_hat(:,n);               end;
    if loop_shock == 6 NXAds_hat   = shocks.NXAds_hat(:,n);            end;
        
    
    ww      = wwM_hat(1:(N-1));
    M_hat   = (reshape(wwM_hat(N:end),[N J]));
    wA      = (1-ww.*sum(bl.Rds0(N,:)))./ sum(bl.Rds0(1,:));        
    w_hat   = [wA;ww];
    
    M_oNNJY = repmat(permute(M_hat,[1 3 2]),[1 N 1]);
    ham     = repmat(w_hat,[1 N J]).^(-pm.theta_nnj) .* Gamma_hat;                                                   
    diff1   = repmat(w_hat,[1 J]) - squeeze(sum( (bl.zeta0.* ham .* repmat(permute(beta_hat,[3 1 2]),[N 1 1])) .* repmat(((w_hat.*sum(bl.Rds0,2)-NXd_hat.*bl.NXd0) ./(bl.Rd0-bl.NXd0))',[N 1 J]) ./ repmat(sum( bl.lambda0.*M_oNNJY.*ham ,1),[N 1 1]),2));
    diff1   = reshape(diff1,[N*J,1]);

    pwrA       = (pm.sigma_nj-1).*(1+pm.theta_nj)./(pm.sigma_nj.*pm.theta_nj);    
    pwrB       = (pm.theta_nj-(pm.sigma_nj-1).*(1-pm.alpha_nj))./(pm.theta_nj.*pm.sigma_nj);
    pwrC       = (pm.sigma_nj-1).*(pm.theta_nj-pm.alpha_nj+1)./(pm.sigma_nj.*pm.theta_nj);
    bl.NXd02D  = repmat(bl.NXd0  ,[1 J]);
    NXAds02D   = repmat(bl.NXAds0 ,[1 J]);
    NXd_hat2D  = repmat(NXd_hat  ,[1 J]);
    NXAds_hat2D= repmat(NXAds_hat,[1 J]);
    eta        = squeeze(sum(-(pwrB-1).*bl.beta0.*bl.NXd02D,2))  - bl.NXAds0;
    eta_p      =   (1./w_hat).*squeeze(sum(-(pwrB-1).*beta_hat.*bl.beta0.*(NXd_hat2D.*bl.NXd02D),2)) - (1./w_hat).*bl.NXAds0.*NXAds_hat;
    psi = squeeze((1-sum(pwrB.*bl.beta0,2)) ./ (1-sum(pwrB.*bl.beta0.*beta_hat,2)));
    diff2      = 1 - psi.* (squeeze(sum(M_hat.*bl.Rds0.*pwrC,2))+eta_p) ./(squeeze(sum(    bl.Rds0.*pwrC,2))+eta  );
                 
    diff = [diff2(2:end);diff1];
    
q = diff;