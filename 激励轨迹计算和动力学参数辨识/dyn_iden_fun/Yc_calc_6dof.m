function Yc = Yc_calc_6dof(dq)
% Friction regressor for 6 joints

dq_min = 0.001 * ones(1,6);
Kc = zeros(1,6);
Kv = zeros(1,6);

for i = 1:6
    if abs(dq(i)) >= dq_min(i)
        Kc(i) = sign(dq(i));
        Kv(i) = dq(i);
    end
end

Yc = zeros(6,12);
for i = 1:6
    col = (i-1)*2 + 1;
    Yc(i, col) = Kc(i);
    Yc(i, col+1) = Kv(i);
end

end
