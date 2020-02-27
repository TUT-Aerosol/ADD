function [out] = dp_davidson(stab_class,x)
% Calculates dispersion parameters (sig_i) following modified Pasquill-Gifford formulae
% (Davidson, J. Air & Waste Man. Assoc., 1990) Here it is assumed that sig_x= sig_y. 
% % Note that the results are NOT valid at distances <100 m from the source.
%
% INPUT: 
% * x - distance from the source (m)
% * stab_class - atmospheric stability class
% 
% OUTPUT:
% * sig_i (i= x, y, z) - dispersion parameters (in m) 
% * dsig_i (i= x, y, z) - derivates of dispersion parameters (in m/m)

out.stab_flag= true; % see below
xx= x/1000.0; % switch to km

switch lower(stab_class)
    case('a')
        y= [209.6	0.8804 -0.006902];
        z= [417.9	2.058	0.2499];
    case('b') 
        y= [154.7	0.8932 -0.006271]; 
        z= [109.8	1.064	0.01163]; 
    case('c') 
        y= [103.3	0.9112 -0.004845]; 
        z= [61.14	0.9147	0]; 
    case('d') 
        y= [68.28	0.9112 -0.004845]; 
        z= [30.38	0.7309	-0.032]; 
    case('e') 
        y= [51.05	0.9112 -0.004845]; 
        z= [21.14	0.6802	-0.04522]; 
    case('f')
        y= [33.96	0.9112 -0.004845]; 
        z= [13.72	0.6584	-0.05367]; 
    otherwise 
        out.stab_flag= false;  % unknown stability class
        disp('** Unknown stability class **'); 
end

if (out.stab_flag == true && all(xx>0))
    out.sig_z= z(1).*(xx.^(z(2)+z(3).*log(xx))); % sig is in m
    out.sig_y= y(1).*(xx.^(y(2)+y(3).*log(xx))); 
    out.sig_x= out.sig_y; 
    out.dsig_z= (z(2)+2.*z(3)).*out.sig_z./x; % dsig is in m / m 
    out.dsig_y= (y(2)+2.*y(3)).*out.sig_y./x; 
    out.dsig_x= out.dsig_y;
else
    out.sig_z= 0; out.sig_y= 0; out.sig_x= 0;
    out.dsig_x= 0.0; out.dsig_y= 0.0; out.dsig_z= 0.0;     
end

end

    
