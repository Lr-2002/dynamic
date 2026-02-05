function [w, dw ,dv ,dvc] = motion_para_clc(R_inv,dq,ddq,w_pre,dw_pre,dv_pre,Po,Pc,axis)
    w = R_inv*w_pre+dq*axis;
    dw = R_inv*dw_pre+cross(R_inv*w_pre,dq*axis)+ddq*axis;
    dv = R_inv*(cross(dw_pre,Po)+cross(w_pre,cross(w_pre,Po))+dv_pre);
    dvc = cross(dw,Pc)+cross(w,cross(w,Pc))+dv;
end
