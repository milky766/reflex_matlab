%% 
data = readtable('20240818_12.xlsx', 'ReadVariableNames', true, 'VariableNamingRule', 'preserve');

%%
% モデルのパラメータ agonist
m = -60.1;   
h_const = 170.1; % 'h'はMATLABで予約語なので別名に
a0 = -16.084; %a_0_newのこと（voltageに対応）
a1 = -4.52271;
a2 = -123.616;
a3 = 3.545457;

%%
% サンプリング時間
dt = 0.01; % 例：0.01秒

% ノイズの標準偏差（センサー仕様や実験データに基づいて設定）
sigma_w_l = 0.01;  % 長さのプロセスノイズ標準偏差（例）
sigma_w_v = 0.1;   % 速度のプロセスノイズ標準偏差（例）
sigma_v_V = 0.05;  % 電圧の観測ノイズ標準偏差（例）

% 共分散行列
Q = [sigma_w_l^2, 0; 0, sigma_w_v^2];
R = sigma_v_V^2;

% データの読み込み
t = data.count*0.01;       % 時間
p = data.pressure_left;       % 圧力
V = data.tension_left;       % 電圧


% データの長さ
N = length(t);

% 状態推定値の初期化
x_est = zeros(2, N);
x_est(:,1) = [0; 0]; % 初期値を適切に設定（例：0, 0）

% 誤差共分散行列の初期化
P = [1, 0; 0, 1]; % 初期誤差共分散（適切に設定）

% フィルタリングループ
for k = 1:N-1
    % サンプリング時間間隔の取得
    Delta_t = t(k+1) - t(k);
    
    % 状態遷移行列の設定
    F = [1, Delta_t; 0, 1];
    
    % 予測ステップ
    x_pred = F * x_est(:,k);
    P_pred = F * P * F' + Q;
    

    % ヤコビアンの計算
    dh_dl = derivative_of_h(x_pred(1), p(k), m, h_const, a0, a1, a2, a3);
    H = [dh_dl, 0];
   
    % 変形量の予測
    d_pred = x_pred(1) - (m * p(k) + h_const);
    
    % 観測関数の予測
    h_pred = (a3 * p(k) * d_pred + a2 * p(k) + a1 * d_pred + a0) * d_pred;
    

    % カルマンゲインの計算
    K = P_pred * H' / (H * P_pred * H' + R);
    
    % 観測残差の計算
    y_k = V(k) - h_pred;
    
    % 状態推定値の更新
    x_est(:,k+1) = x_pred + K * y_k;
    
    % 誤差共分散行列の更新
    P = (eye(2) - K * H) * P_pred;
end

%%
% 推定結果の抽出
l_est = x_est(1, :); % 推定された長さ
v_est = x_est(2, :); % 推定された速度

% ms_left_modelデータの取得
ms_left_model = data.ms_left_model;

% ms_speed_left_modelデータの取得
ms_speed_left_model = data.ms_speed_left_model;


% 長さの推定値と ms_left_model の比較プロット
figure;
plot(t, l_est, 'r-', 'LineWidth', 2); % 推定された長さのプロット
hold on;
plot(t, ms_left_model, 'g--', 'LineWidth', 2); % ms_left_modelのプロット
xlabel('時間 [s]');
ylabel('長さ [単位]');
title('推定された長さと ms\_left\_model の比較');
legend('推定された長さ', 'ms\_left\_model');
grid on;

% 速度の推定値と ms_speed_left_model の比較プロット
figure;
plot(t, v_est, 'b-', 'LineWidth', 2); % 推定された速度のプロット
hold on;
plot(t, ms_speed_left_model, 'm--', 'LineWidth', 2); % ms_speed_left_modelのプロット
xlabel('時間 [s]');
ylabel('速度 [単位/s]');
title('推定された速度と ms\_speed\_left\_model の比較');
legend('推定された速度', 'ms\_speed\_left\_model');
grid on;
