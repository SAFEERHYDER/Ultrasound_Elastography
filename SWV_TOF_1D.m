
%   Given the local displacement shear waves data
%   Returns the estimated  local shear wave velocity (SWV) at each pixel
%   location


            %    Author  : SAFEER HYDER
            %    Updated : 21-04-2016


%         parameter.w                     : Lateral Kernal
%         parameter.Filter                : 1 = Yes, 0 = No
%         parameter.Filter_coefficients   : Filter coefficients
%         parameter.Filter_Order          : Filter Order
%         parameter.Interpolation_order   : Interpolation factor
%         parameter.taper_factor          : window tapering factor
%         parameter.pixel_in_m            : Lateral pixel in m
%         parameter.PRP                   : Frame period



function [SWVEL, C, x_axis]  = SWV_TOF_1D(DATA, parameter)

clear C SWVEL;
[depth, lat, ~] = size(DATA);
lat_start   = 1;
lat_end     = lat - parameter.W;
depth_start = 1;
depth_end   = depth;

% SWVEL  = zeros(depth, radial - parameter.W);


if (~parameter.point && parameter.slice)
    depth_start = parameter.depth;
    depth_end   = parameter.depth;
end

if (parameter.point && parameter.slice)
    depth_start  = parameter.depth;
    depth_end    = parameter.depth;
    lat_start    = parameter.lat;
    lat_end      = parameter.lat;
end


for axial = depth_start : depth_end
    for lateral = lat_start: lat_end
        
        
        
        S1 = -squeeze( DATA(axial, lateral, :))*parameter.ps*1E6;
        S2 = -squeeze( DATA(axial, lateral + parameter.W, :))*parameter.ps*1E6;
        
        if parameter.Mean
            S1 = S1 - mean(S1);
            S2 = S2 - mean(S2);
        end
        
        if parameter.Filter
            S1 = conv(S1, parameter.Filter_coefficients);
            S1 = S1(parameter.Filter_Order/2 +1:end-parameter.Filter_Order/2);
            S2 = conv(S2, parameter.Filter_coefficients);
            S2 = S2(parameter.Filter_Order/2 +1:end-parameter.Filter_Order/2);
        end
        
        
        if parameter.Interp
            S1 = Interp_disp(S1, parameter.Interpolation_order);
            S2 = Interp_disp(S2, parameter.Interpolation_order);
        end
        
        
        if parameter.window
            window = tukeywin(numel(S1), parameter.taper_factor  );
            S1 = S1.*window';
            S2 = S2.*window';
        end
        
        
        [coeff, Delay_Lag_Interp] = xcorr( S1,  S2, 'coeff');
        [C(axial, lateral), C_index ] = max(coeff);
        Delay_Lag_True = Delay_Lag_Interp(C_index) / parameter.Interpolation_order;
        
        distance =  parameter.pixel_in_m*parameter.W;
        time     =  abs(parameter.PRP*Delay_Lag_True);
        SWVEL (axial, lateral) = distance / time;
        
        tp = (parameter.PRP/parameter.Interpolation_order)*1E3;
        t = tp:tp:numel(S1)*tp;
        
        if parameter.plot
            figure();
            subplot(2,1,1);
            plot(t, S1, 'blue');
            xlabel('\bf slow time (ms)');
            ylabel('\bf displacement (\mu m)');
            hold on; plot(t, S2, 'red');
            legend(['Lateral = ', num2str(lateral*parameter.pixel_in_m*1E3), ' mm'],['Lateral = ',num2str((lateral+parameter.W)*parameter.pixel_in_m*1E3), ' mm']);
            title(['Axial = ', num2str(axial*0.5), ' mm'])
            subplot(2,1,2);
            plot(Delay_Lag_Interp*tp, coeff); xlim([-1.5 +1.5]);
            xlabel('\bf lag (ms)');
            ylabel('\bf correlation coefficient');
        end
    end
end




x_axis = - (lateral/2)*parameter.pixel_in_m:parameter.pixel_in_m:(lateral/2)*parameter.pixel_in_m;
x_axis = x_axis(ceil(parameter.W/2) : end- floor(parameter.W/2) );
x_axis = x_axis*1E3;

if (parameter.point && parameter.slice)
    x_axis = parameter.lat;
end

end




