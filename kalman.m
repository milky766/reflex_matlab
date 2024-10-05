% モデルのパラメータ
%agonist
m = -60.1;    % 実験で決定された値
h_const = 170.1; % 'h'はMATLABで予約語なので別名に
a0 = 0.91102;
a1 = ...;
a2 = ...;
a3 = ...;

% サンプリング時間
dt = 0.01; % 例：0.01秒

% プロセスノイズの標準偏差
sigma_w_l = ...;
sigma_w_v = ...;

% 観測ノイズの標準偏差
sigma_v_f = ...;

% プロセスノイズ共分散行列
Q = [sigma_w_l^2, 0; 0, sigma_w_v^2];

% 観測ノイズ共分散
R = sigma_v_f^2;
