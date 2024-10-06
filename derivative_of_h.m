function dh_dl = derivative_of_h(l, p, m, h_const, a0, a1, a2, a3)
    % 変形量の計算
    d = l - (m * p + h_const);
    
    % 観測関数のdに関する偏微分
    dh_dd = 2*(a3 * p + a1) * d + (a2 * p + a0);
    
    % dのlに関する偏微分
    dd_dl = 1;
    
    % hのlに関する偏微分
    dh_dl = dh_dd * dd_dl;
end
