function Flux=LaxFriedFlux(UL,UR,k,g,h)
%compute eigenvalues using formula nablau f((UL+UR)/2)
U=(UL+UR)*0.5;
lamda1=U(2,:)./U(1,:) + sqrt(g*(U(1,:)));
lamda2=U(2,:)./U(1,:) - sqrt(g*(U(1,:)));
lam=max((abs(lamda1)),(abs(lamda2)));
%Definition of the flux
Flux=0.5*([UL(2,:);UL(2,:).^2./UL(1,:) + 0.5*g*UL(1,:).^2]+[UR(2,:);UR(2,:).^2./UR(1,:) + 0.5*g*UR(1,:).^2])-(1/((2)))*(UR-UL).*lam;
end