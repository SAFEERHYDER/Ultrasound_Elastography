clc
clear all
close all
% path (path,'E:\MECE\UltraSound\TEST\cyst_phantom_test\Field_II_combined' )

field_init;

%% Create a transducer 

fs=100e6; % Sampling frequency [Hz]
Ts=1/fs;
c=1540; % Speed of sound [m/s]
height=6e-3; % Height of element [m]
L=c*Ts/2; % L is the distance between each 2 points in the z direction
jjj=1;

% for f0=10e6:4e6:30e6
    
f0=4e6;

lambda=c/f0; 
pitch=0.3048e-3;

N_elements=96;

kerf= pitch/10; % Distance between transducer elements [m]
width=pitch-kerf; % Width of element [m]

%%%%%%%%%%%%%%Initial electronic focus%%%%%%%%%%%%%%%
focus=[0 0 100]; 
Th = xdc_linear_array (N_elements, width, height, kerf, 1,1, focus);
Th2 = xdc_linear_array (N_elements, width, height, kerf, 1,1, focus);

BW=2*f0;

impulse_response=sin(2*pi*f0*(0:1/fs:2/BW));
impulse_response=impulse_response.*hanning(max(size(impulse_response)))';


xdc_impulse (Th, impulse_response);
xdc_impulse (Th2, impulse_response);

BW=2e6/f0;

%  for BW= 0.1 : 0.1 : 2
    tc=gauspuls('cutoff' ,f0 ,BW,[],-60); % the BW here specifies only the duration
    ... of the total pulse but not the amount of energy inside the envelope, i.e. it pads zeros to reach this duration but no change on the signal.
    tt=-tc:Ts:tc;
    excitation=gauspuls(tt, f0, BW); % the BW here specifies the duration of the pulse or the number of sycles in the pulse.

xdc_excitation (Th, excitation);
Half_signal=(length(excitation) + 2*length(impulse_response)-2)*0.5;
    
%%%%%%%%%%%%%%%%%%%% Scattering Points phantom created manually %%%%%%%%%%%%%%%%%%
% [positions, amp] = cyst_ph(100000);
positions=[0 0 5]/1000;
amp=1;
imaging_start = 1/1000; % (A-3)/1000;
imaging_depth = 10/1000; % 6/1000;
image_width = 10/1000;
dx=lambda/20;

yy=zeros(25000,10000); % initial matrix for addition
W=(N_elements*pitch); % -kerf
L=c*Ts/2; % L is the distance between each 2 points in the z direction
z_start = imaging_start;
z_stop = imaging_start+imaging_depth;
N_activeT=N_elements;
Tapo=hamming(N_activeT)';
xdc_apodization (Th, 0, Tapo);
n=1:N_elements;


    %%%%%%%%%%%%%%%%%% SETTING THE DELAYS %%%%%%%%%%%%%%%%
   
    N_angles=1;
    
    for thetad=0:3:0 
        delays=n*pitch*sind(thetad)/c;
        delays=delays-min(delays);
        xdc_focus_times(Th, 0, delays);

        %%%%%%%%% using calc_scat %%%%%%%%%%%%%%%%
        v=zeros(8000,N_activeT);
        
        for i=1:N_activeT
            Rapo=[zeros(1, i-1) 1 zeros(1, N_elements-i)];
            xdc_apodization (Th2, 0, Rapo);
            [v1,t1]=calc_scat(Th, Th2, positions, amp);
            v1_zeropadded = [zeros(1,round(t1*fs))'; v1];
            v(1:length(v1_zeropadded),i) = v1_zeropadded;
        end

time_array=1:size(v,1);
x=-image_width/2 : dx : image_width/2 ;
z=z_start: L : z_stop ;
ff= zeros(length(z), length(x));
[X, Z]=meshgrid(x,z);
d1= Z; % *cosd(theta_d)+X*sind(theta_d)+(W/2)*sind(abs(theta_d)); % the distance from the transmitter to the points    
elements_counter=1;

for xj= -(N_elements/2)+0.5 : (N_elements/2)-0.5  % loop for the receiving elements
    RF_address= (( d1+sqrt(Z.^2+(X-xj*pitch).^2) ) ./ (c*Ts))+Half_signal;
    ff= ff+ interp1(time_array,v(:,elements_counter),RF_address, 'linear', 0);
    elements_counter=elements_counter+1;
end

%%%%%%%%%%%%%% Hilbert Transform %%%%%%%%%%%%%%
        yy(1:size(ff,1),1:size(ff,2))=yy(1:size(ff,1),1:size(ff,2))+ff;
  
    end
    
    ab1=yy(1:size(ff,1),1:size(ff,2))/N_angles;
    env=abs(hilbert(ab1/max(max(ab1))));
    
    %%%%%%%%%%%%%%%%%%% Display %%%%%%%%%%%%%%%%%%%%%%
    
    env_dB=20*log10(env);
    
    depth=z_start:L:z_stop;
    x=-image_width/2:dx:image_width/2;
    figure;
    imagesc(x*1000, depth*1000, env_dB, [-50 0])
    xlabel('Lateral distance [mm]')
    ylabel('Depth [mm]')
%     title('Four scattering points at 5, 10, 15 and 20cm depths'); %(['Aperture Width= ' num2str(apwidth*1e2) 'cm'])
    axis('image')
    colormap(gray(128));
    colorbar;
    
% %     axis square
% %%
%     [M N]=size(env_dB);  % env_dB instead of J (J and env_dB have the same size)
%     Sz=(M/2); % careful zainab!! the size is important
%     Sx=(N/6);
%     for a=1:M
%         for b=1:N
%             if (env_dB(a,b)<-50)
%                 env_dB(a,b)=-50;
%             elseif (env_dB(a,b)>0)
%                 env_dB(a,b)=0;
%             end
%         end
%     end
%     
%     mm=1;
%     b=env_dB(round(mm*Sz),:)
%     
%     figure;
%     plot((1:N)/Sx-3,env_dB(round(mm*Sz),:) );   % for the point at 30mm (mm is 5 as the image starts from 25mm).
%     [MLW3dB FWHM FWHDR PSL DML] = calculate_performance_metrics(env_dB(round(mm*Sz),:), dx, 50); 
% %     figure;
% %     plot((1:M)/Sz+27,env_dB(:, round(N/2)) );
% %     xlabel('Lateral distance [mm]')
% %     ylabel('Intensity [dB]')
% %     title(['the Lateral resolution at the 30 mm depth, N='])
%     LatR_6(jjj)=FWHM;
%     jjj=jjj+1;
% % end
% 
% AAA=1:0.5:30;
% figure; plot(AAA,LatR_6*1e3);
%     grid on
%     xlabel('Frequency [MHz]');
%     ylabel('Lateral Resolution [mm]');
%     return
%     %%%%%%%%%%%%%%%% To calculate Lateral Resolution:
%     aaa=round(mm*Sz)+16;
%     [MLW3dB FWHM FWHDR PSL DML] = calculate_performance_metrics(env_dB(aaa,:), dx, 50); 
%     
% %     LatR_3(jjj)=MLW3dB; 
% 
% % %%%%%%%%%%%%%%%%%%% To caulcuate Axial Resolution:
% %      [MLW3dB FWHM FWHDR PSL DML] = calculate_performance_metrics(env_dB(:, round(N/2)) , L , 50); 
% % %     
% % % %     AxR_3(jjj)=MLW3dB; 
% %     AxR_6(jjj)=FWHM;
%     jjj=jjj+1;
% % %     xxx(jjj)=positions(3)/(f0^2 * N_elements^2*apwidth
% % end
% 
% return;
% % % % %% the values of latR vs. lambda, when attenuation exists:
% % % % f=1e6:1e6:10e6;
% % % % c=1540;
% % % % lmda=c./f;
% % % % LR_att=[3.1644  1.64  1.0452  .7494  .5797  .4481  .3786  .3163  .2826  .251]*1e-3; % this is in mm.
% % % % % the values of latR vs. lambda, when attenuation does not exist;
% % % % LR_no_att=[0.00301734564388659,0.00149305140083842,0.000916865948677397,0.000649257554424107,0.000497127666657826,0.000373464528193191,0.000318201456750865,0.000272340612938370,0.000244478488144671,0.000218848021590983];% this is in mm.
% % % % figure; plot(f/1e6, LR_att/1e-3,f/1e6, LR_att/1e-3,'s' , f/1e6, LR_no_att/1e-3, f/1e6, LR_no_att/1e-3,'s')
% % % % legend('with attenuation','','no attenuation',''); grid on
% % % % %%%%%%%%%%%%%%%%%%%%%%%%%
% % % % %%
% % % % apwidth=10e-3:5e-3:60e-3;
% % % % LR_Att=[0.00154700207846632,0.00103149795062346,0.000795714473470410,0.000660832128928334,0.000574229650101793,0.000510694205980678,0.000470306002363360,0.000438325381342294,0.000414176198808855,0.000395229803006250,0.000381495714335108];
% % % % LR_noAtt=[0.00135568841074934,0.000903713208894717,0.000694469358591331,0.000574105950180470,0.000497127666657826,0.000438071813021517,0.000402244429105654,0.000375616007320863,0.000355543041057476,0.000341400101994423,0.000329315950037724];
% % % % figure; plot(apwidth/1e-3, LR_Att/1e-3, apwidth/1e-3, LR_noAtt)
% % % % legend('with attenuation','no attenuation'); grid on
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
