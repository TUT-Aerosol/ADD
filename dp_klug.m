function [out] = dp_klug(stab_class, x)
% Calculates dispersion parameters (sig_i) following Klug(1969) formulae
% (Seinfeld, 1985, Table 14.2, second to last column). Here it assumed that sig_x= sig_y. 
% Note that the results are NOT valid at distances <100 m from the source.
%
% INPUT: 
% * x - distance from the source (m)
% * stab_class - atmospheric stability class
% 
% OUTPUT:
% * sig_i (i= x,y,z) - dispersion parameters (in m)
% * dsig_i (i= x, y, z) - derivates of dispersion parameters (in m/m)

out.stab_flag= true;

switch lower(stab_class)
    case('a')
        z= [0.017  1.380];
        y= [0.469  0.903];
    case('b') 
        z= [0.072 1.021]; 
        y= [0.306 0.885]; 
    case('c') 
        z= [0.076 0.879]; 
        y= [0.230 0.855]; 
    case('d') 
        z= [0.140 0.727]; 
        y= [0.219 0.764]; 
    case('e') 
        z= [0.217 0.610]; 
        y= [0.237 0.691]; 
    case('f')
        z= [0.262 0.500]; 
        y= [0.273 0.594]; 
    otherwise 
        out.stab_flag= false;  % unknown stability class
%        disp('** Unknown stability class **') 
end

if (out.stab_flag == true && x > 0)
    out.sig_z= z(1)*(x^z(2));
    out.sig_y= y(1)*(x^y(2));
    out.sig_x= out.sig_y; 
    out.dsig_z= z(2)*out.sig_z/x; 
    out.dsig_y= y(2)*out.sig_y/x;
    out.dsig_x= out.dsig_y;
else
    out.stab_flag= false; 
    out.sig_z= 0; out.sig_y= 0; out.sig_x= 0;
    out.dsig_z= 0.0; out.dsig_y= 0.0; out.dsig_x= 0.0;
end

end

    
