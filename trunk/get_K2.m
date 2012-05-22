function K = get_K2( BW, T )
%GET_K2 Get coeffs for 2nd-order locking loop

K = [0; 0];

K(2) = 2*16/9*BW^2; % For continious model
K(1) = sqrt(2*K(2));
K = K*T; % For discrete model

end

