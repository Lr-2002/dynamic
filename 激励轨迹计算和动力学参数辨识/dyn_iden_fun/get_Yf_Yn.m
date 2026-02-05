function [Yf, Yn] = get_Yf_Yn(id, R, H, A, Yf_next, Yn_next, Po)
c1=(id-1)*10+1;c2=(id-1)*10+10;
Yf=zeros(3,80);Yf(:,c1:c2)=H;
Yf=Yf+R*Yf_next;
Yn=zeros(3,80);Yn(:,c1:c2)=A;
Yn=Yn+R*Yn_next+getS(Po)*R*Yf_next;

end
