

% parameter.slice                         = 15;
% parameter.point                         = 0;
% parameter.plot                          = 1;
% [swvel_L, swvelcof_L, x_axis]            = SWV_TOF_1D_Main(axdisp_MID_L, parameter);
% 




function [SWVEL, C, x_axis] = SWV_TOF_1D_Main(DATA3D, slice, point, parameter)

% close all;

% velocity for 2D Grid
if ~(parameter.slice  || parameter.point)
    [SWVEL, C, x_axis]  = SWV_TOF_1D(DATA3D, parameter);
    
    
% velocity for single slice
elseif (~parameter.point && parameter.slice)
    parameter.depth = slice;
    [SWVEL, C, x_axis]  = SWV_TOF_1D(DATA3D, parameter);
    
    
% velocity at single point
elseif (parameter.slice && parameter.point)
    parameter.depth = slice;
    parameter.lat   = point;
    [SWVEL, C, x_axis]  = SWV_TOF_1D(DATA3D, parameter);
end



end


