function [X, Y, Z] = get_vector_XYZ( r, lambda, i, u )
%GET_VECTOR_XYZ Get vector XYZ by his 3 angles and length

xyz = U3(-lambda)*U1(-i)*U3(-u)*[r; 0; 0];

X = xyz(1);
Y = xyz(2);
Z = xyz(3);

end

