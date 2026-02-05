function R = RotZ(theta)

    ct = cos(theta);
    st = sin(theta);

    R = [
        ct  -st  0
        st   ct  0
        0    0   1
        ];

end

