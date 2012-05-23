function K = get_K3( BW, T )
%GET_K3 Get coeffs for 3rd-order locking loop

K = [0; 0; 0];

K(3) = (1.2*BW)^3; % For continious model
K(2) = 2*(K(3))^(2/3);
K(1) = 2*(K(3))^(1/3);

K = K*T; % For discrete model

end

