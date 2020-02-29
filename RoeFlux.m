function Flux=RoeFlux(UL,UR,g)
%compute eigenvalues using formula nablau f((UL+UR)/2)

zl=UL./sqrt(UL(1,:));
zr=UR./sqrt(UR(1,:));
x=((zl(2,:)+zr(2,:))*0.5)./((zl(1,:)+zr(1,:))*0.5)+ sqrt((0.5*g*(zl(1,:).^2 + zr(1,:).^2)));
y=((zl(2,:)+zr(2,:))*0.5)./((zl(1,:)+zr(1,:))*0.5)- sqrt((0.5*g*(zl(1,:).^2 + zr(1,:).^2)));
Flux=zeros(2,size(UR,2));
for i=1:size(UR,2)
    U=[UR(1,i)-UL(1,i);UR(2,i)-UL(2,i)];
    S=[1 1; x(i) y(i)];
    Sinv=[y(i)/(y(i)-x(i)) -(y(i)-x(i))^-1; -x(i)/(y(i)-x(i))   y(i)^-1+x(i)/(y(i)*(y(i)-x(i)))];
    Lambda=[x(i) 0;0 y(i)];
    absA=S*abs(Lambda)*Sinv;
    %Aplus=S*[max(x(i),0) 0;0 max(y(i),0)]*Sinv;
    %Aminus=S*[max(-x(i),0) 0;0 max(-y(i),0)]*Sinv;
    Flux(:,i)=0.5*([UL(2,i);UL(2,i)^2/UL(1,i) + 0.5*g*UL(1,i)^2]+[UR(2,i);UR(2,i)^2/UR(1,i) + 0.5*g*UR(1,i)^2])-(1/((2)))*(absA)*(U);
end

end