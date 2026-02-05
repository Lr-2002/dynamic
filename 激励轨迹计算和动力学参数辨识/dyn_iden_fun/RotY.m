function R = RotY(theta)

    ct = cos(theta);
    st = sin(theta);

    R = [
        ct  0   st
        0   1   0
       -st  0   ct
       ];

end
