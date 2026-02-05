function I = getInertiaMatrix(Ic_vec)

    Ixx = Ic_vec(1); Iyy = Ic_vec(2); Izz = Ic_vec(3);
    Iyz = Ic_vec(4); Ixz = Ic_vec(5); Ixy = Ic_vec(6);

    I = [Ixx, Ixy, Ixz;
        Ixy, Iyy, Iyz;
        Ixz, Iyz, Izz];

end
