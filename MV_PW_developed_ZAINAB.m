%%% func_MV_PW_developed.m
%   Action: 
%       Performs time domain Minimum Variance beamforming (MV) on the RF data of a
%       plane wave steered with (theta_d) degree.
%
%   Usage:
%       [Beamformed_data] = func_MV_PW_developed (RF_data, theta_d, z_start, z_stop,
%       image_width, W, N_elements, pitch, c, Ts, td, Lp, dx)

%   Inputs:
%       RF_data is the received RF data.
%       theta_d is the steering angle of the plane wave in degree, can be
%       a positive or a negative value.
%       z_start is the imaging depth starting point (in meters).
%       z_stop is the imaging depth ending point (in meters).
%       image_width is the required width of the produced image.
%       W is the total width of the transducer (in meters), where: W=N_elements x pitch.
%       lambda is the wavelength.
%       N_elements is the total number of elements in the transducer that
%       are used to transmit the ultrasound beam.
%       pitch is the distance between the centres of two adjacent elements
%       (in meter).
%       c is the sound speed im m/sec.
%       Ts is the sampling time.
%       td is the length of segment taken from each element during the
%       beamforming. The centre of this segment resembles the response
%       received from the focal point. To retain axial resolution, td
%       should not exceed the pulse length (which is the convolution
%       between the excitation signal and the 2-way impulse response).
%       (td is always an odd number).
%       Lp is the number of elements in each subarray. Increasing Lp gives
%       higher resolution on the cost of lowering robustness. Decreasing Lp
%       reduces the resolution and improves the robustness, where Lp=1
%       resembles DAS beamformer with uniform aperture shading (unity apodization).
%       dx is the resolution (or step) in x direction.
%
%   Returns:
%           Beamformed_data: The frame resulted from MV beamforming operation.
%
%   Notes:
%      In MV beamforming, a vector of length (td) is selected from each
%      element's RF-data, where the response from the focal point is at the
%      centre point of this vector.
%      For each point (x,z), the delays applied to the RF data are
%      calculated by the time required for the signal to travel from the
%      transmitter to the field point and back to the receiving element as follows:
%      t= t_transmit + t_receive
%      t_transmit= [ z.cos(theta_d)+x.sin(theta_d)+0.5.W.sin(|theta_d|) ]/c
%      t_receive(for element j)= sqrt(z^2 + (x-xj)^2)/c
%      (where xj is the axial distance of the j-th element)

%      After performing the delays, adpative weight is calculated for each
%      point according to the equation that ensures minimizing the output
%      power while preserving a unity gain in the focal point. This equation is: 
%      w=(inv(R).e)/(conj(e)'.inv(R).e)
%      where e is a vector of ones. R is the covariance matrix.
%      After calculating the weight, the output value is calculated using:
%      B=conj(w)'.summation(Gp)/P
%      where P is the number of subarrays and Gp is the subarray.

%       Diagonal Loading can be used to increase the robutness of the
%       beamformer by adding a constant value to the diagonal of the
%       covariance matrix using the following equation:
%       R=R+ delta.trace{R}.I
%       where delta is 1/(DL.Lp), DL is the diagonal loading coefficient,
%       I is the identity matrix.

%      The resulted image dimensions are specified by the user by z_start,
%      z_stop and image_width. For compounding 3 angles for example,
%      this function will be called 3 times and the three reults are
%      averaged, normalized to the maximum value, envelope detected
%      by hilbert transform and then converted to dB.

% References:
% 1) J-F Synnevag et al., "Adaptive beamforming applied to medical
% ultrasound imaging", Ultrasonics, Ferroelectrics and Frequency Control,
% IEEE Transactions on, 54(8):1606-1613, 2007.
% 2) I. K. Holfort et al., "Broadband minimum variance beamforming for
% ultrasound imaging", Ultrasonics, Ferroelectrics and Frequency Control, 
% IEEE Transactions on, 56(2):314-325, 2009.

%       To decrease the processing time (when only Plane-Wave imaging is
%       being performed), comment line (104) and uncomment line (105) for d1.



function [Beamformed_data]= func_MV_PW_developed(RF_data, theta_d, z_start , z_stop  , image_width, dx , W, N_elements, pitch , c , Ts , td , Lp)
tic
L=c*Ts/2;
P=N_elements-Lp+1;  % the number of subarrays
e=ones(Lp,1);  
tdd= (td-1)/2;
time_array=1:size(RF_data,1);
Z=z_start: L :z_stop;
X=-image_width/2: dx: image_width/2;
tdd_matrix=repmat((-tdd:tdd) ,length(Z),1); % size of this matrix is (length(Z) x td)

%%%%% Memory allocation:
Ym=zeros(length(Z), td, N_elements);
Beamformed_data=zeros(length(Z), length(X));

adrs_x=1;
for x=-image_width/2:dx:image_width/2
    d1= Z.*cosd(theta_d)+x*sind(theta_d)+(W/2)*sind(abs(theta_d));% the distance from the transmitter to the points that have the same lateral distance of x.
%     d1=Z; % d1 becomes only Z when transmitting with angle 0. 
    xjj=1;
    for xj=-(N_elements/2)+0.5:(N_elements/2)-0.5
        address_vector=(d1+sqrt(Z.^2+(x-(xj*pitch))^2))/c/Ts;
        address_matrix= repmat(address_vector',1,td);             % size of this matrix is (length(Z) x td)
        final_matrix=address_matrix + tdd_matrix;                 % a matrix containing a td-length address vector for each point of x lateral distance.
        Ym(:,:,xjj)=interp1(time_array, RF_data(:,xjj), final_matrix, 'linear', 0); % Ym contains a td-length vector for each point of x lateral distance for each element. 
        xjj=xjj+1;
    end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %%%%%%%%%%% now calculating the weight (w) for each point:
    for adrs_z=1:length(Z)
        R=zeros(Lp,Lp); % initializing covariance matrix
        sum_Gp=zeros(Lp,td);                    
        for loopG = 1 : P
            Gp=(shiftdim(Ym(adrs_z , : , loopG:loopG+Lp-1 ),1))'; % Converting from 3d to 2d.
            R=R+(Gp * conj(Gp'));
            sum_Gp=sum_Gp+Gp;
        end 
        R=R/P;
%         R=R+(1/(10*Lp)) *trace(R)*eye(Lp); %applying Diagonal Loading to R before calculating the weight.
        w=(R \ e) / (conj(e') * (R \ e)) ; % calculating the weight.
        B= conj(w') * (sum_Gp / P);
        Beamformed_data(adrs_z,adrs_x)=B((td+1)/2); % taking the central point of the output vector.
    end
    adrs_x=adrs_x+1;
end
toc
