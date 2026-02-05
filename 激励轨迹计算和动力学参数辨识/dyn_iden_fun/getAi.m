function Ai=getAi(wi,dwi,dvi)
Ai(1:3,1:6)=getK(dwi)+getS(wi)*getK(wi);
Ai(1:3,7:9)=-getS(dvi);
Ai(1:3,10)=zeros(3,1);
end
