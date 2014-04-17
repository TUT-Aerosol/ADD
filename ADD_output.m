function [] = ADD_output(out, dpout)
%ADD_output.m - plots results for agglo_disp_driv 

% choose scaling for the 1st plot: 
% axis_scaling= 'linear'; 
axis_scaling= 'semilog';

switch(axis_scaling)
    case('semilog')
        figure(01);
        subplot(1,3,1); semilogy(out.tc_ts, out.Ntot_ts,'LineWidth',2);
        grid on; axis tight;
        xlabel('Time (s)'); ylabel('Number concentration (/cm^3)');
        % title('Number concentration');
        subplot(1,3,2); semilogy(out.tc_ts, out.phi_ts,'LineWidth',2);
        grid on; axis tight; 
        xlabel('Time (s)'); ylabel('Total aerosol volume concentration (m^3/cm^3)');
        % title('Volume concentration');
        subplot(1,3,3); plot(out.tc_ts, out.va_ts/out.v_a0,'LineWidth',2);
        grid on; axis tight; 
        xlabel('Time (s)'); ylabel('Relative particle volume (V/V_0)');
        % title('Relative particle volume');
    case('linear')
        figure(01);
        subplot(1,3,1); plot(out.tc_ts, out.Ntot_ts,'LineWidth',2);
        grid on; axis tight;
        xlabel('Time (s)'); ylabel('Number concentration (/cm^3)');
        % title('Number concentration');
        subplot(1,3,2); semilogy(out.tc_ts, out.phi_ts,'LineWidth',2);
        grid on; axis tight;
        xlabel('Time (s)'); ylabel('Total aerosol volume concentration (m^3/cm^3)');
        % title('Volume concentration');
        subplot(1,3,3); plot(out.tc_ts, out.va_ts/out.v_a0,'LineWidth',2);
        grid on; axis tight;
        xlabel('Time (s)'); ylabel('Relative particle volume (V/V_0)');
        % title('Relative particle volume');
end % end switch
 
figure(02);  
subplot(1,2,1); plot(out.x_ts, out.sigx_ts, 'LineWidth',2);
grid on; axis tight;
xlabel('x (m)'); ylabel('sigma_x (m)');
subplot(1,2,2); plot(out.x_ts, out.sigz_ts, 'LineWidth',2);
grid on; axis tight; 
xlabel('x (m)'); ylabel('sigma_z (m)');

figure(03); 
plot(out.tc_ts, out.x_ts, 'r', out.tc_ts, out.y_ts, 'g', out.tc_ts, out.z_ts, 'b', 'LineWidth',2);
% title('Coordinates');
grid on; 
xlabel('Time (s)'); ylabel('Distance (m)');
legend('x', 'y', 'z');

figure(04); 
semilogy(out.tc_ts, out.Kl_ts, 'LineWidth', 2);
grid on; axis tight;
xlabel('Time (s)'); ylabel('Loss coefficient (/s)');

figure(05); 
subplot(2,3,1); plot(out.tc_ts, out.cx_ts, 'LineWidth',2);
grid on;  axis tight;
xlabel('Time (s)'); ylabel('c_x');
subplot(2,3,2); plot(out.tc_ts, out.cy_ts, 'LineWidth',2);
grid on;  axis tight;
xlabel('Time (s)'); ylabel('c_y');
subplot(2,3,3); plot(out.tc_ts, out.cz_ts, 'LineWidth',2);
grid on;  axis tight;
xlabel('Time (s)'); ylabel('c_z');
subplot(2,3,4); plot(out.tc_ts, out.dcx_ts, 'LineWidth',2);
grid on;  axis tight;
xlabel('Time (s)'); ylabel('dc_x/dt (/s)');
subplot(2,3,5); plot(out.tc_ts, out.dcy_ts, 'LineWidth',2);
grid on;  axis tight;
xlabel('Time (s))'); ylabel('dc_y/dt (/s)');
subplot(2,3,6); plot(out.tc_ts, out.dcz_ts, 'LineWidth',2);
grid on;  axis tight;
xlabel('Time (s)'); ylabel('dc_z/dt (/s)');


% for testing: 
% figure(06); 
% semilogy(out.x_ts, out.Ctest_ts, 'LineWidth', 2);
% grid on; axis tight;
% xlabel('Distance (m)'); ylabel('Q');

end

