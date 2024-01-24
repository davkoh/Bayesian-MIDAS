function [beta,hyper_params,data] = sample_beta(data,hyper_params)


%% Unpack Data

% Data
yhat = data.Y;
X = data.X;
tX = X';
iOh = data.iOh;
gl_param_expand_diag_inv = data.gl_param_expand_diag_inv;
K = data.K;
tau_shape_const = data.tau_shape_const;
stable_const = data.stable_const;
G = data.G;
grp_size_cs = data.grp_size_cs;
grp_size = data.grp_size;
grp_idx = data.grp_idx;
a = data.a;
b = data.b;


% Hyperparameters
tau_sq = hyper_params.tau_sq;
nu = hyper_params.nu;
gamma_sq = hyper_params.gamma_sq;
lambda_sq = hyper_params.lambda_sq;


for gg  = 1:K
    gl_param_expand_diag_inv(gg) = 1.0 / (tau_sq * gamma_sq(grp_idx(gg)) * lambda_sq(gg));
end

    beta_tmp = tX*iOh*X + sparse(diag(gl_param_expand_diag_inv)) + 1e-10;
    
    beta = tX*iOh*yhat + chol(beta_tmp,'lower')*randn(K,1);%beta = (1.0 / sigma_sq) * tX * (Y) + chol(beta_tmp,'lower')*randn(M,1);
   
    beta = beta_tmp\beta;


%% Draw tau^2
tau_rate_const = sum(beta.^2.*gl_param_expand_diag_inv);
tau_sq = 1.0 / gamrnd(tau_shape_const, 1.0 / (tau_sq * tau_rate_const / 2.0 + 1.0 / nu));

%% Draw gamma_g^2/lambda^2_gj
for j = 1:G

    % Sample gamma_g^2
    stable_psi = 0;
    if j == 1
    start_tmp = 1;
    end_tmp = grp_size_cs(j);
    else
    start_tmp = grp_size_cs(j-1)+1;
    end_tmp = grp_size_cs(j);
    end

    stable_psi = sum(beta(start_tmp:end_tmp).^2./lambda_sq(start_tmp:end_tmp));   
    stable_psi = stable_psi./tau_sq;
    stable_psi = max(stable_psi,stable_const);
    gamma_sq(j) = 1/gigrnd(grp_size(j)/2-a(j),stable_psi, 2, 1); %%%% Watch out for the p variable here.

    

    % Sample lambda^2_gj
    for i = 1:grp_size(j)
        lambda_sq(start_tmp+i-1) = 1.0 / gamrnd(b(j) + 0.5,...
            1.0 / (1 + (beta(start_tmp + i-1)^2) / (2.0 * tau_sq * gamma_sq(j))));
    end
end

%% Draw nu
nu = 1.0 / gamrnd(1, 1 / ((1 / tau_sq) ));


data.gl_param_expand_diag_inv = gl_param_expand_diag_inv;
hyper_params.tau_sq = tau_sq;
hyper_params.nu = nu;
hyper_params.gamma_sq = gamma_sq;
hyper_params.lambda_sq = lambda_sq;

