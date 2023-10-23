function output_tf = tf_from_ss(A,B,C,D)
%tf_from_ss: Makes a transfer function from A, B, C, and D matrices found
%in the state-space representation
    s = tf('s');
    phi = (s*eye(size(A,1))-A)\1;
    output_tf = C*phi*B+D;
end