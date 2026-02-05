function Hi=getHi(wi,dwi,dvi)
Hi(1:3,7:9)=getS(dwi)+getS(wi)*getS(wi);
Hi(1:3,10)=dvi;
Hi(1:3,1:6)=zeros(3,6);
end
