clc; clear all; close all;

shot = 88567; t1 = 0.3000; t2 = 0.3012;

fmagin = get_scd_mems(shot,'fmagin');
fmagin = fmagin.fmagin;
bm_rt = fmagin.bm_rt.Data; ff_rt = fmagin.ff_rt.Data; uf_rt = fmagin.uf_rt.Data;
t = fmagin.t.Time;

t_mask = (t >= t1) & (t <= t2);
t = t(t_mask);
bm_rt = bm_rt(t_mask, :);
ff_rt = ff_rt(t_mask, :);
uf_rt = uf_rt(t_mask, :);

[L,LX] = liuqe(shot, t);
bm_lx = LX.bm_rt'; ff_lx = LX.ff_rt'; uf_lx = LX.uf_rt';
e_bm = abs(bm_rt - bm_lx); e_ff = abs(ff_rt - ff_lx); e_uf = abs(uf_rt - uf_lx);

for i = 1:length(t)
    fprintf('Time: %.6f \n', t(i));
    fprintf('Err bm_rt: min:%.6f, max:%.6f, mean:%.6f \n', min(e_bm(i,:)), max(e_bm(i,:)), mean(e_bm(i,:)));
    fprintf('Err ff_rt: min:%.6f, max:%.6f, mean:%.6f \n', min(e_ff(i,:)), max(e_ff(i,:)), mean(e_ff(i,:)));
    fprintf('Err uf_rt: min:%.6f, max:%.6f, mean:%.6f \n', min(e_uf(i,:)), max(e_uf(i,:)), mean(e_uf(i,:)));   
end