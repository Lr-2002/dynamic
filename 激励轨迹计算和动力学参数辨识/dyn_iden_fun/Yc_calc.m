function Yc = Yc_calc(dq)
% 有关摩擦力参数辨识的矩阵

dq1=dq(1); dq2=dq(2); dq3=dq(3);dq4=dq(4); dq5=dq(5); dq6=dq(6); dq7=dq(7); dq8=dq(8);
dq_min_1 = 0.001; dq_min_2 = 0.001; dq_min_3 = 0.001; dq_min_8 = 0.001;
dq_min_4 = 0.001; dq_min_5 = 0.001; dq_min_6 = 0.001; dq_min_7 = 0.001;

Kc1 = 0; Kv1 = 0; Kc2 = 0; Kv2 = 0; Kc3 = 0; Kv3 = 0;
Kc4 = 0; Kv4 = 0; Kc5 = 0; Kv5 = 0; Kc6 = 0; Kv6 = 0;

if abs(dq1) >= dq_min_1
    Kc1 = sign(dq1);
    Kv1 = dq1;
else
    Kc1 = 0;
    Kv1 = 0;
end

if abs(dq2) >= dq_min_2
    Kc2 = sign(dq2);
    Kv2 = dq2;
else
    Kc2 = 0;
    Kv2 = 0;
end

if abs(dq3) >= dq_min_3
    Kc3 = sign(dq3);
    Kv3 = dq3;
else
    Kc3 = 0;
    Kv3 = 0;
end

if abs(dq4) >= dq_min_4
    Kc4 = sign(dq4);
    Kv4 = dq4;
else
    Kc4 = 0;
    Kv4 = 0;
end

if abs(dq5) >= dq_min_5
    Kc5 = sign(dq5);
    Kv5 = dq5;
else
    Kc5 = 0;
    Kv5 = 0;
end

if abs(dq6) >= dq_min_6
    Kc6 = sign(dq6);
    Kv6 = dq6;
else
    Kc6 = 0;
    Kv6 = 0;
end

if abs(dq7) >= dq_min_7
    Kc7 = sign(dq7);
    Kv7 = dq7;
else
    Kc7 = 0;
    Kv7 = 0;
end

if abs(dq8) >= dq_min_8
    Kc8 = sign(dq8);
    Kv8 = dq8;
else
    Kc8 = 0;
    Kv8 = 0;
end

Yc = [Kc1 Kv1 0 0 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 Kc2 Kv2 0 0 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 Kc3 Kv3 0 0 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 Kc4 Kv4 0 0 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 Kc5 Kv5 0 0 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 Kc6 Kv6 0 0 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 Kc7 Kv7 0 0;
      0 0 0 0 0 0 0 0 0 0 0 0 0 0 Kc8 Kv8];

end
