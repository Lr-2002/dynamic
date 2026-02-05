function R = RotX(theta)

    ct = cos(theta);
    st = sin(theta);

    R = [
        1   0    0
        0   ct  -st
        0   st   ct
        ];

end
