
    %       CIRS Homogeneous Phantom Study


  
    
%%     Directory of the Ultrasound DATA and Elastography Matlab Functions  


close all;
clear; 
cd('D:\UARPII Data\UARPII\UARP - DPB - Experiments\28_06_2016_CPWI_SWEI_EXPERIMENT');
addpath(genpath('D:\UARPII Data\UARPII\UARP - DPB - Experiments\28_06_2016_CPWI_SWEI_EXPERIMENT'));


%%  Load shear wave displacement data


load('AP1-DF', 'z_axis', 'x_axis', 'axdisp_MID', 'BF', 'z', 'ps');
Exp1.Name = 'a) DPB-AP1';
Exp1.Axdisp = axdisp_MID(:,:,1:130);


load('AP2-DF', 'z_axis', 'x_axis', 'axdisp_MID', 'BF', 'z', 'ps');
Exp2.Name = 'b) DPB-AP2';
Exp2.Axdisp = axdisp_MID(:,:,1:130);


load('AP3-DF', 'z_axis', 'x_axis', 'axdisp_MID', 'BF', 'z', 'ps');
Exp3.Name = 'c) DPB-AP3';
Exp3.Axdisp = axdisp_MID(:,:,1:130);


load('AP4-DF', 'z_axis', 'x_axis', 'axdisp_MID', 'BF', 'z', 'ps');
Exp4.Name = 'd) DPB-AP4';
Exp4.Axdisp = axdisp_MID(:,:,1:130);


load('SSI-LEFT', 'z_axis', 'x_axis', 'axdisp_MID', 'BF', 'z', 'ps');
Exp5.Name = 'e) SSI-LEFT';
Exp5.Axdisp = axdisp_MID(:,:,1:130);


load('SSI-RIGHT', 'z_axis', 'x_axis', 'axdisp_MID', 'BF', 'z', 'ps');
Exp6.Name = 'f) SSI-RIGHT';
Exp6.Axdisp = axdisp_MID(:,:,1:130);



%% Parameters


% Low Pass Filter
set.Interp_f                        = 10;
set.PRF                             = 10E3;
set.show_plots                      = 0;
set.l_cutoff                        = 10;
set.u_cutoff                        = 500;
set.wl                              = set.l_cutoff / (set.PRF/2);
set.wu                              = set.u_cutoff / (set.PRF/2);
set.order                           = 10;
set.Filter_coefficients             = fir1(set.order, [set.wl set.wu]);

% Arrival-time algorithm parameters
set.PRF                             = set.PRF;
set.do_window                       = 1;
set.do_lowpass                      = 1;
set.do_demean                       = 1;
set.interp_f                        = set.Interp_f;
set.lateral_step                    = BF.lateral_step;
set.dlat                            = 8;


% Image depth dimension
z_axis                              = Interp_disp(z, 2);
axial_step                          = (z_axis(2) - z_axis(1))*1E-3 ;
set.z_stop                          = 40e-3;
z_start_index                       = 15;      % 3 mm
z_stop_index                        = ceil(set.z_stop / axial_step) - 5;


%% Event#01 : TWO-BEAM-Ap1-SF

set.z_stop                          = 40e-3;
z_stop_index                        = ceil(set.z_stop / axial_step) - 5;



% Interp in the axial direction 2 x 
for frame = 1:size(Exp1.Axdisp, 3)
    for lateral = 1:size(Exp1.Axdisp, 2)
        DISP_DATA_1(:, lateral, frame) = Interp_disp(Exp1.Axdisp(:, lateral, frame), 2);
    end
end 


% spatial filtering : removing outliers
for frame = 1:size(DISP_DATA_1, 3)
    DISP_DATA_2(:,:,frame) = medfilt2(DISP_DATA_1(:, :, frame), [3 3]);
end



% LR shear wave 
axdisp_MID_LR  = direct_filter_3D(DISP_DATA_2,     1, 0.05);
[swvel_R, swvel_R_coeff]   = Est_local_swv_inv_ORG(axdisp_MID_LR,   [], [], set);
Exp1.B1        = abs(swvel_R(:, set.dlat +1:end - set.dlat));


% figure(1);     
% % subplot(3,4,subfig); 
% % imagesc(x_axis(set.dlat +1:end - set.dlat), z_axis(z_start_index:z_stop_index) ,(fliplr(Exp1.B1(1:end, 1:end).^2)*1.050*3), [1 100]);
% imagesc(Exp1.B1, [2 6]);
% axis image
% colormap jet
% colorbar; c = colorbar;
% ylabel (c,'\bf Youngs modulus (kPa)', 'FontSize', 9);
% xlabel('\bf Lateral (mm)', 'Fontsize', 9);
% ylabel('\bf Depth   (mm)', 'FontSize', 9);
% title([num2str(Exp1.Name)]);


% RL shear waves 
axdisp_MID_RL  = direct_filter_3D(DISP_DATA_2,     0, 0.05);
[swvel_L, swvel_L_coeff]   = Est_local_swv_inv_ORG(axdisp_MID_RL,   [], [], set);
Exp1.B2        = abs(swvel_L(:, set.dlat +1:end - set.dlat));


% figure(2); 
% % subplot(3,4,subfig); 
% % imagesc(x_axis(set.dlat +1:end - set.dlat), z_axis(z_start_index:z_stop_index) ,(fliplr(Exp1.B2(1:end, 1:end).^2)*1.050*3), [1 100]);
% imagesc(Exp1.B2, [2 6]);
% axis image
% colormap jet
% colorbar; c = colorbar;
% ylabel (c,'\bf Youngs modulus (kPa)', 'FontSize', 9);
% xlabel('\bf Lateral (mm)', 'Fontsize', 9);
% ylabel('\bf Depth   (mm)', 'FontSize', 9);
% title([num2str(Exp1.Name)]);


% LR-RL shear waves 
RL_fused_image =  [ Exp1.B2(:, 1:121) Exp1.B1(:, 122:end) ];
RL_coeff_fused1 = [swvel_R_coeff(:, 1:122) swvel_L_coeff(:, 123:end)];
RL_coeff_fused2 = RL_coeff_fused1(:, set.dlat +1:end - set.dlat);

% for axial = 1:size(RL_coeff_fused2,1)
%     for lateral = 1:size(RL_coeff_fused2, 2)
%         if RL_coeff_fused2(axial, lateral) >= 0.6
%             Mask_coef (axial, lateral) = 1;
%         else 
%             Mask_coef (axial, lateral) = 0;
%         end 
%     end 
% end 

% RL_filt_masked = RL_fused_image.*Mask_coef; 
Exp1.RL_fused_image_filt = fliplr(medfilt2(RL_fused_image, [4 4]));
figure(3); 
subplot(4,2,1); 
imagesc(x_axis(set.dlat +1:end - set.dlat), z_axis(z_start_index:z_stop_index), ((Exp1.RL_fused_image_filt(z_start_index:z_stop_index, :).^2)*1.050*3), [1 60]);
% imagesc((Exp1.RL_fused_image_filt.^2)*1.050*3, [1 100]);
axis image
colormap jet
colorbar; c = colorbar;
ylabel (c,'\bf Youngs modulus (kPa)', 'FontSize', 9);
xlabel('\bf Lateral (mm)', 'Fontsize', 9);
ylabel('\bf Depth (mm)', 'FontSize', 9);
title([num2str(Exp1.Name)]);


% ROI statistics 
shear_mod_map   = (Exp1.RL_fused_image_filt(z_start_index:z_stop_index, :).^2)*1.050*3;
counter1 = 1;
for axial = 1:size(shear_mod_map,1)
    for lateral = 1:size(shear_mod_map,2)
        if shear_mod_map(axial, lateral) >= 50
            shear_mod_map(axial, lateral) =  22.5;
            counter1 = counter1 + 1;
        end 
    end 
end 

figure(3); 
subplot(4,2,1); 
imagesc(x_axis(set.dlat +1:end - set.dlat), z_axis(z_start_index:z_stop_index), shear_mod_map, [1 60]);
% imagesc((Exp1.RL_fused_image_filt.^2)*1.050*3, [1 100]);
axis image
colormap jet
colorbar; c = colorbar;
ylabel (c,'\bf Youngs modulus (kPa)', 'FontSize', 9);
xlabel('\bf Lateral (mm)', 'Fontsize', 9);
ylabel('\bf Depth (mm)', 'FontSize', 9);
title([num2str(Exp1.Name)]);


Exp1.Mean       = mean(mean(shear_mod_map)); 
Exp1.STD        = std(std(shear_mod_map));


%% Event#02 : TWO-BEAM-Ap1-DF

set.z_stop                          = 40e-3;
z_stop_index                        = ceil(set.z_stop / axial_step) - 5;

subfig = 2 ;


% Interp in the axial direction 2 x 
for frame = 1:size(Exp2.Axdisp, 3)
    for lateral = 1:size(Exp2.Axdisp, 2)
        DISP_DATA_1(:, lateral, frame) = Interp_disp(Exp2.Axdisp(:, lateral, frame), 2);
    end
end 


% spatial filtering : removing outliers
for frame = 1:size(DISP_DATA_1, 3)
    DISP_DATA_2(:,:,frame) = medfilt2(DISP_DATA_1(:, :, frame), [3 3]);
end


% LR shear wave 
axdisp_MID_LR  = direct_filter_3D(DISP_DATA_2,     1, 0.05);
[swvel_R, swvel_R_coeff]   = Est_local_swv_inv_ORG(axdisp_MID_LR,   [], [], set);
Exp2.B1        = abs(swvel_R(:, set.dlat +1:end - set.dlat));


% figure(1);     
% % subplot(3,4,subfig); 
% % imagesc(x_axis(set.dlat +1:end - set.dlat), z_axis(z_start_index:z_stop_index) ,(fliplr(Exp2.B1(1:end, 1:end).^2)*1.050*3), [1 100]);
% imagesc(Exp2.B1, [2 6]);
% axis image
% colormap jet
% colorbar; c = colorbar;
% ylabel (c,'\bf Youngs modulus (kPa)', 'FontSize', 9);
% xlabel('\bf Lateral (mm)', 'Fontsize', 9);
% ylabel('\bf Depth   (mm)', 'FontSize', 9);
% title([num2str(Exp2.Name)]);


% RL shear waves 
axdisp_MID_RL  = direct_filter_3D(DISP_DATA_2,     0, 0.05);
[swvel_L, swvel_L_coeff]   = Est_local_swv_inv_ORG(axdisp_MID_RL,   [], [], set);
Exp2.B2        = abs(swvel_L(:, set.dlat +1:end - set.dlat));


% figure(2); 
% % subplot(3,4,subfig); 
% % imagesc(x_axis(set.dlat +1:end - set.dlat), z_axis(z_start_index:z_stop_index) ,(fliplr(Exp2.B2(1:end, 1:end).^2)*1.050*3), [1 100]);
% imagesc(Exp2.B2, [2 6]);
% axis image
% colormap jet
% colorbar; c = colorbar;
% ylabel (c,'\bf Youngs modulus (kPa)', 'FontSize', 9);
% xlabel('\bf Lateral (mm)', 'Fontsize', 9);
% ylabel('\bf Depth   (mm)', 'FontSize', 9);
% title([num2str(Exp2.Name)]);


% % LR-RL shear waves 
RL_fused_image =  [ Exp2.B2(:, 1:121) Exp2.B1(:, 122:end) ];
% RL_coeff_fused1 = [swvel_R_coeff(:, 1:122) swvel_L_coeff(:, 123:end)];
% RL_coeff_fused2 = RL_coeff_fused1(:, set.dlat +1:end - set.dlat);
% 
% for axial = 1:size(RL_coeff_fused2,1)
%     for lateral = 1:size(RL_coeff_fused2, 2)
%         if RL_coeff_fused2(axial, lateral) >= 0.6
%             Mask_coef (axial, lateral) = 1;
%         else 
%             Mask_coef (axial, lateral) = 0;
%         end 
%     end 
% end 


% RL_filt_masked = RL_fused_image.*Mask_coef; 
Exp2.RL_fused_image_filt = fliplr(medfilt2(RL_fused_image, [4 4]));



% ROI statistics 
shear_mod_map   = (Exp2.RL_fused_image_filt(z_start_index:z_stop_index, :).^2)*1.050*3;
counter2 = 1;
for axial = 1:size(shear_mod_map,1)
    for lateral = 1:size(shear_mod_map,2)
        if shear_mod_map(axial, lateral) >= 50
            shear_mod_map(axial, lateral) =  22.5;
            counter2 = counter2 + 1;
        end 
    end 
end


figure(3); 
subplot(4,2,2); 
imagesc(x_axis(set.dlat +1:end - set.dlat), z_axis(z_start_index:z_stop_index), shear_mod_map, [1 60]);
% imagesc((Exp1.RL_fused_image_filt.^2)*1.050*3, [1 100]);
axis image
colormap jet
colorbar; c = colorbar;
ylabel (c,'\bf Youngs modulus (kPa)', 'FontSize', 9);
xlabel('\bf Lateral (mm)', 'Fontsize', 9);
ylabel('\bf Depth (mm)', 'FontSize', 9);
title([num2str(Exp2.Name)]);

Exp2.Mean       = mean(mean(shear_mod_map)); 
Exp2.STD        = std(std(shear_mod_map));


%% Event#03 : TWO-BEAM-Ap2-SF


set.z_stop                          = 40e-3;
z_stop_index                        = ceil(set.z_stop / axial_step) - 5;


% Interp in the axial direction 2 x 
for frame = 1:size(Exp3.Axdisp, 3)
    for lateral = 1:size(Exp3.Axdisp, 2)
        DISP_DATA_1(:, lateral, frame) = Interp_disp(Exp3.Axdisp(:, lateral, frame), 2);
    end
end 


% spatial filtering : removing outliers
for frame = 1:size(DISP_DATA_1, 3)
    DISP_DATA_2(:,:,frame) = medfilt2(DISP_DATA_1(:, :, frame), [3 3]);
end


% LR shear wave 
axdisp_MID_LR  = direct_filter_3D(DISP_DATA_2,     1, 0.05);
[swvel_R, swvel_R_coeff]   = Est_local_swv_inv_ORG(axdisp_MID_LR,   [], [], set);
Exp3.B1        = abs(swvel_R(:, set.dlat +1:end - set.dlat));


% figure(1);     
% subplot(3,1,1); 
% imagesc(x_axis(set.dlat +1:end - set.dlat), z_axis(z_start_index:z_stop_index) ,(fliplr(Exp3.B1(1:end, 1:end).^2)*1.050*3), [1 100]);
% colormap jet
% colorbar; c = colorbar;
% ylabel (c,'\bf Youngs modulus (kPa)', 'FontSize', 9);
% xlabel('\bf Lateral (mm)', 'Fontsize', 9);
% ylabel('\bf Depth   (mm)', 'FontSize', 9);
% axis image


% RL shear waves 
axdisp_MID_RL  = direct_filter_3D(DISP_DATA_2,     0, 0.05);
[swvel_L, swvel_L_coeff]   = Est_local_swv_inv_ORG(axdisp_MID_RL,   [], [], set);
Exp3.B2        = abs(swvel_L(:, set.dlat +1:end - set.dlat));


% figure(1); 
% subplot(3,1,2); 
% imagesc(x_axis(set.dlat +1:end - set.dlat), z_axis(z_start_index:z_stop_index) ,(fliplr(Exp3.B2(1:end, 1:end).^2)*1.050*3), [1 100]);
% axis image
% colormap jet
% colorbar; c = colorbar;
% ylabel (c,'\bf Youngs modulus (kPa)', 'FontSize', 9);
% xlabel('\bf Lateral (mm)', 'Fontsize', 9);
% ylabel('\bf Depth   (mm)', 'FontSize', 9);


% % LR-RL shear waves 
RL_fused_image =  [ Exp3.B2(:, 1:121)  Exp3.B1(:, 122:end) ];
% RL_coeff_fused1 = [swvel_R_coeff(:, 1:122) swvel_L_coeff(:, 123:end)];
% RL_coeff_fused2 = RL_coeff_fused1(:, set.dlat +1:end - set.dlat);

% for axial = 1:size(RL_coeff_fused2,1)
%     for lateral = 1:size(RL_coeff_fused2, 2)
%         if RL_coeff_fused2(axial, lateral) >= 0.6
%             Mask_coef (axial, lateral) = 1;
%         else 
%             Mask_coef (axial, lateral) = 0;
%         end 
%     end 
% end 


% RL_filt_masked = RL_fused_image.*Mask_coef; 
Exp3.RL_fused_image_filt = fliplr(medfilt2(RL_fused_image, [4 4]));




% ROI statistics 
shear_mod_map   = (Exp3.RL_fused_image_filt(z_start_index:z_stop_index, :).^2)*1.050*3;
counter3 = 1;
for axial = 1:size(shear_mod_map,1)
    for lateral = 1:size(shear_mod_map,2)
        if shear_mod_map(axial, lateral) >= 50
            shear_mod_map(axial, lateral) =  22.5;
            counter3 = counter3 + 1;
        end 
    end 
end

figure(3); 
subplot(4,2,3); 
imagesc(x_axis(set.dlat +1:end - set.dlat), z_axis(z_start_index:z_stop_index), shear_mod_map, [1 60]);
% imagesc((Exp1.RL_fused_image_filt.^2)*1.050*3, [1 100]);
axis image
colormap jet
colorbar; c = colorbar;
ylabel (c,'\bf Youngs modulus (kPa)', 'FontSize', 9);
xlabel('\bf Lateral (mm)', 'Fontsize', 9);
ylabel('\bf Depth (mm)', 'FontSize', 9);
title([num2str(Exp3.Name)]);



Exp3.Mean       = mean(mean(shear_mod_map));
Exp3.STD        = std(std(shear_mod_map));


%% Event#04 : TWO-BEAM-Ap4-DF

set.z_stop                          = 40e-3;
z_stop_index                        = ceil(set.z_stop / axial_step) - 5;

subfig = 4;

% Interp in the axial direction 2 x 
for frame = 1:size(Exp4.Axdisp, 3)
    for lateral = 1:size(Exp4.Axdisp, 2)
        DISP_DATA_1(:, lateral, frame) = Interp_disp(Exp4.Axdisp(:, lateral, frame), 2);
    end
end 


% spatial filtering : removing outliers
for frame = 1:size(DISP_DATA_1, 3)
    DISP_DATA_2(:,:,frame) = medfilt2(DISP_DATA_1(:, :, frame), [3 3]);
end


% LR shear wave 
axdisp_MID_LR  = direct_filter_3D(DISP_DATA_2,     1, 0.05);
[swvel_R, swvel_R_coeff]   = Est_local_swv_inv_ORG(axdisp_MID_LR,   [], [], set);
Exp4.B1        = abs(swvel_R(:, set.dlat +1:end - set.dlat));


% figure(1);     
% % subplot(3,4,subfig); 
% % imagesc(x_axis(set.dlat +1:end - set.dlat), z_axis(z_start_index:z_stop_index) ,(fliplr(Exp4.B1(1:end, 1:end).^2)*1.050*3), [1 100]);
% imagesc(Exp4.B1, [2 6]);
% axis image
% colormap jet
% colorbar; c = colorbar;
% ylabel (c,'\bf Youngs modulus (kPa)', 'FontSize', 9);
% xlabel('\bf Lateral (mm)', 'Fontsize', 9);
% ylabel('\bf Depth   (mm)', 'FontSize', 9);
% title([num2str(Exp4.Name)]);


% RL shear waves 
axdisp_MID_RL  = direct_filter_3D(DISP_DATA_2,     0, 0.05);
[swvel_L, swvel_L_coeff]   = Est_local_swv_inv_ORG(axdisp_MID_RL,   [], [], set);
Exp4.B2        = abs(swvel_L(:, set.dlat +1:end - set.dlat));


% figure(2); 
% % subplot(3,4,subfig); 
% imagesc(x_axis(set.dlat +1:end - set.dlat), z_axis(z_start_index:z_stop_index) ,(fliplr(Exp4.B2(1:end, 1:end).^2)*1.050*3), [1 100]);
% imagesc(Exp4.B2, [2 6]);
% axis image
% colormap jet
% colorbar; c = colorbar;
% ylabel (c,'\bf Youngs modulus (kPa)', 'FontSize', 9);
% xlabel('\bf Lateral (mm)', 'Fontsize', 9);
% ylabel('\bf Depth   (mm)', 'FontSize', 9);
% title([num2str(Exp4.Name)]);


% LR-RL shear waves 
RL_fused_image =  [ Exp4.B2(:, 1:121)  Exp4.B1(:, 122:end) ];
Exp4.RL_fused_image_filt = fliplr(medfilt2(RL_fused_image, [4 4]));



% ROI statistics 
shear_mod_map   = (Exp4.RL_fused_image_filt(z_start_index:z_stop_index, :).^2)*1.050*3 ;
counter4 = 1;
for axial = 1:size(shear_mod_map,1)
    for lateral = 1:size(shear_mod_map,2)
        if shear_mod_map(axial, lateral) >= 50
            shear_mod_map(axial, lateral) =  22.5;
            counter4 = counter4 + 1;
        end 
    end 
end

figure(3); 
subplot(4,2,4); 
imagesc(x_axis(set.dlat +1:end - set.dlat), z_axis(z_start_index:z_stop_index), shear_mod_map, [1 60]);
% imagesc((Exp4.RL_fused_image_filt.^2)*1.050*3, [1 100]);
axis image
colormap jet
colorbar; c = colorbar;
ylabel (c,'\bf Youngs modulus (kPa)', 'FontSize', 9);
xlabel('\bf Lateral (mm)', 'Fontsize', 9);
ylabel('\bf Depth (mm)', 'FontSize', 9);
title([num2str(Exp4.Name)]);


Exp4.Mean       = mean(mean(shear_mod_map));
Exp4.STD        = std(std(shear_mod_map));




%% Event#05 : SSI

set.z_stop                          = 40e-3;
z_stop_index                        = ceil(set.z_stop / axial_step) - 5;


% Interp in the axial direction 2 x 
for frame = 1:size(Exp5.Axdisp, 3)
    for lateral = 1:size(Exp5.Axdisp, 2)
        DISP_DATA_1(:, lateral, frame) = Interp_disp(Exp5.Axdisp(:, lateral, frame), 2);
    end
end 


% spatial filtering : removing outliers
for frame = 1:size(DISP_DATA_1, 3)
    DISP_DATA_2(:,:,frame) = medfilt2(DISP_DATA_1(:, :, frame), [3 3]);
end


% LR shear wave
axdisp_MID_RL  = direct_filter_3D(DISP_DATA_2,     1, 0.05);
[swvel_R, swvel_R_coeff]   = Est_local_swv_inv_ORG(axdisp_MID_RL,   [], [], set);
Exp5.B1        = abs(swvel_R(:, set.dlat +1:end - set.dlat));
Exp5.B1filt = medfilt2(Exp5.B1, [4,4]);


% 
figure(3);     
subplot(4,2,5); 
imagesc(x_axis(set.dlat +1:end - set.dlat), z_axis(z_start_index:z_stop_index) ,fliplr((Exp5.B1filt(z_start_index:z_stop_index, :).^2)*1.050*3), [1 60]);
axis image
colormap jet
colorbar; c = colorbar;
ylabel (c,'\bf Youngs modulus (kPa)', 'FontSize', 9);
xlabel('\bf Lateral (mm)', 'Fontsize', 9);
ylabel('\bf Depth   (mm)', 'FontSize', 9);
title([num2str(Exp5.Name)]);

%% SSI


% Interp in the axial direction 2 x 
for frame = 1:size(Exp6.Axdisp, 3)
    for lateral = 1:size(Exp6.Axdisp, 2)
        DISP_DATA_1(:, lateral, frame) = Interp_disp(Exp6.Axdisp(:, lateral, frame), 2);
    end
end 


% Spatial filtering : removing outliers
for frame = 1:size(DISP_DATA_1, 3)
    DISP_DATA_2(:,:,frame) = medfilt2(DISP_DATA_1(:, :, frame), [3 3]);
end


% LR shear wave
axdisp_MID_RL  = direct_filter_3D(DISP_DATA_2,     0, 0.05);
[swvel_R, swvel_R_coeff]   = Est_local_swv_inv_ORG(axdisp_MID_RL,   [], [], set);
Exp6.B1        = abs(swvel_R(:, set.dlat +1:end - set.dlat));
Exp6.B1filt = medfilt2(Exp6.B1, [4 4]);


figure(3);     
subplot(4,2,6); 
imagesc(x_axis(set.dlat +1:end - set.dlat), z_axis(z_start_index:z_stop_index) ,fliplr((Exp6.B1filt(z_start_index:z_stop_index, :).^2)*1.050*3), [1 60]);
axis image
colormap jet
colorbar; c = colorbar;
ylabel (c,'\bf Youngs modulus (kPa)', 'FontSize', 9);
xlabel('\bf Lateral (mm)', 'Fontsize', 9);
ylabel('\bf Depth   (mm)', 'FontSize', 9);
title([num2str(Exp6.Name)]);



%% SSI combined


% LR-RL shear waves 
RL_fused_image =  [ Exp6.B1(:, 1:121)  Exp5.B1(:, 122:end) ];
Exp56.RL_fused_image_filt = fliplr(medfilt2(RL_fused_image, [4 4]));


% ROI statistics 
shear_mod_map   = (Exp56.RL_fused_image_filt(z_start_index:z_stop_index, :).^2)*1.050*3 ;
counter56 = 1;
for axial = 1:size(shear_mod_map,1)
    for lateral = 1:size(shear_mod_map,2)
        if shear_mod_map(axial, lateral) >= 50
            shear_mod_map(axial, lateral) =  22.5;
            counter56 = counter56 + 1;
        end 
    end 
end 



figure(3); 
subplot(4,2,7); 
imagesc(x_axis(set.dlat +1:end - set.dlat), z_axis(z_start_index:z_stop_index), shear_mod_map, [1 60]);
% imagesc((Exp5.RL_fused_image_filt.^2)*1.050*3, [1 100]);
axis image
colormap jet
colorbar; c = colorbar;
ylabel (c,'\bf Youngs modulus (kPa)', 'FontSize', 9);
xlabel('\bf Lateral (mm)', 'Fontsize', 9);
ylabel('\bf Depth (mm)', 'FontSize', 9);
title('g) SSI-Combined');


Exp56.Mean       = mean(mean(shear_mod_map)); 
Exp56.STD        = std(std(shear_mod_map));


%% statistics plots....

Means = [Exp1.Mean Exp2.Mean Exp3.Mean Exp4.Mean Exp56.Mean 22.5]';
STDs  = [Exp1.STD Exp2.STD Exp3.STD Exp4.STD Exp56.STD 5]';

figure(5);
barwitherr(STDs, Means);

  






% %% save Data and Figures
% 
% cd('N:\Faculty-of-Engineering\Research\I3S\Ultrasound\PhD\Safeer\Thesis\SAFEER\Chapter5\images');
% 
% figure(3);
% savefig('DPB-SSI-Youngs-Modulus-Maps-Homo-1');
% 
% 
% 
% 
% h = gcf;
% fig.PaperPositionMode = 'auto';
% print('DPB-SSI-Youngs-Modulus-Maps-Homo-1','-deps','-r400');
% print('DPB-SSI-Youngs-Modulus-Maps-Homo-1','-dpdf','-r400');
% 
% 
% %%
% 
% cd('N:\Faculty-of-Engineering\Research\I3S\Ultrasound\PhD\Safeer\UARPII\UARP - DPB - Experiments\28_06_2016_CPWI_SWEI_EXPERIMENT');
% 
% figure(3);
% savefig('DPB-SSI-Youngs-Modulus-Maps-Homo-1');
% 
% 
% h = gcf;
% fig.PaperPositionMode = 'auto';
% print('DPB-SSI-Youngs-Modulus-Maps-Homo-1','-deps','-r400');
% print('DPB-SSI-Youngs-Modulus-Maps-Homo-1','-dpdf','-r400');





