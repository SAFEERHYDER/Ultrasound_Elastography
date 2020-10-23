% Interpolation function (for quick changing of interpolation method)
% Interplation schemes: Linear, cubic, Spline 

function disp_p = Interp_disp( disp, interp_f )
%disp_p = interp1( disp, 1:1/interp_f:numel(disp), ’linear’ );
%disp_p = interp1( disp, 1:1/interp_f:numel(disp), ’spline’ );
disp_p = interp1( disp, 1:1/interp_f:numel(disp), 'pchip');
%disp_p = interpft( disp, numel(disp)*interp_f );
end
