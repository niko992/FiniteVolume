function U_ext=apply_bc(U,bc)
switch bc
    case 'Periodic'
        U_ext=[U(:,end),U,U(:,1)];
    case 'Open'
        U_ext=[U(:,1),U,U(:,end)];
end
return