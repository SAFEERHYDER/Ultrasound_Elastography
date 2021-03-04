clc
clear all
tic
field_init;
    f0=5e6; % Transducer center frequency [Hz]
    fs=100e6; % Sampling frequency [Hz]
    Ts=1/fs;
    c=1540; %1632.4; %1601.6; % Speed of sound [m/s]
    lambda=c/f0; % Wavelength [m]
    height=6e-3; % Height of element [m]
    pitch=lambda/2;
    kerf=pitch/10; %  Distance between transducer elements [m]
    width=pitch-kerf; % Width of element [m]    
    N_elements=128; % Number of elements

    focus=[0 0 100]; % Initial electronic focus
    Th = xdc_linear_array (N_elements, width, height, kerf, 1,1, focus);
    Th2 = xdc_linear_array (N_elements, width, height, kerf, 1,1, focus);
    BW=1*f0;
    impulse_response=sin(2*pi*f0*(0:Ts:2/BW));
    impulse_response=impulse_response.*hanning(max(size(impulse_response)))';
% % %     tc=gauspuls('cutoff' ,f0 ,1,[],-60); % this is to specify the time duration of the gaussian pulse.
% % %     tt=-tc:Ts:tc;
% % %     impulse_response=gauspuls(tt, f0, 1);
    xdc_impulse (Th, impulse_response);
    xdc_impulse (Th2, impulse_response);
    excitation=sin(2*pi*f0*(0:Ts:2/BW));
    excitation=excitation.*hanning(max(size(excitation)))';
    % excitation=gauspuls((-1/BW:1/fs:1/BW),f0,BW/f0);
% %     tc=gauspuls('cutoff' ,f0 ,0.5,[],-60); this is to specify the time
% duration of the gaussian pulse.
% %     tt=-tc:Ts:tc;
% %     excitation=gauspuls(tt, f0, 0.5);
% % % plot(tt,excitation);
xdc_excitation (Th, excitation);
LengthConv=(length(excitation) + 2*length(impulse_response)-2)*0.5;

dx=lambda/10;
imaging_start = 32/1000;
imaging_depth = 2/1000; % 30/1000; 
image_width=5e-3;
%%
L=c*Ts/2; % L is the distance between each 2 points in the z direction
z_start = imaging_start;
z_stop = imaging_start+imaging_depth;
X=-image_width/2: dx: image_width/2;
Z=z_start: L : z_stop;  
% td_matrix=LengthConv*2;
%%%%%%%%%%%%%% Scattering Points %%%%%%%%%%%%%%%
% [positions, amp, CM] = cyst_phantommm(20000, X, Z);
positions=[0 0 33]/1000; % 0 0 35; 0 0 40; 0 0 45; -3 0 45; 3 0 45; 0 0 50]/1000;
amp=1; %ones(7,1);
% ppp=load('Cphantom'); % Cphantom is the structure that contains both positions and amp, and is stored in the same folder.
    td=LengthConv*2; % this is one pulse length
    tdd=(td-1)/2;
%     Lp_matrix=round(N_elements/4);
%     for Lp_Counter=1:1  % the number of elements in each subarray
        Lp=N_elements/2; %Lp_matrix(Lp_Counter);
        P=N_elements-Lp+1;  % the number of subarrays
        e=ones(Lp,1);  % the steering vector. 
yy=zeros(25000,10000); % initial matrix for addition
W=(N_elements*pitch); % -kerf

        Tapo=hamming(N_elements)';
        xdc_apodization (Th, 0, Tapo);
        
 %%       %%%%%%%%% SETTING THE DELAYS %%%%%%%%%%%
    n=1:N_elements;
    N_angles=1;
    for theta_d=0:5:0
        delays=n*pitch*sind(theta_d)/c;
        delays=delays-min(delays);
        xdc_focus_times(Th, 0, delays);

        [v1,times]=calc_scat_multi(Th, Th2, positions, amp);
        v1_zeropadded= [zeros(round(times*fs) , N_elements); v1];
        v(1:length(v1_zeropadded),:)=v1_zeropadded;
        RF_data=awgn( v, 60, 'measured');
        %%%%%%%%%%% Beamforming %%%%%%%%%%%%
        time_array=1:size(RF_data,1);
        Z=z_start: L :z_stop;
        X=-image_width/2: dx: image_width/2;
        tdd_matrix=repmat((-tdd:tdd) ,length(Z),1); % size of this matrix is (length(Z) x td)
        %%%%% Memory allocation:
        Ym=zeros(length(Z), td, N_elements);
        Beamformed_data=zeros(length(Z), length(X));
%         cf_matrix=zeros(length(Z), length(X));
        adrs_x=1;
        for x=-image_width/2:dx:image_width/2
            d1= Z.*cosd(theta_d)+x*sind(theta_d)+(W/2)*sind(abs(theta_d));% the distance from the transmitter to the points that have the same lateral distance of x.
%             d1=Z; % d1 becomes only Z when transmitting with 0 angle. 
            xjj=1;
            for xj=-(N_elements/2)+0.5:(N_elements/2)-0.5
                address_vector=( (d1+sqrt(Z.^2+(x-(xj*pitch))^2))/c/Ts ) +LengthConv;
                address_matrix= repmat(address_vector',1,td);  % size of this matrix is (length(Z) x td)
                final_matrix=address_matrix + tdd_matrix;% a matrix with a td-length address vector for each point of x lateral distance.
                Ym(:,:,xjj)=interp1(time_array,RF_data(:,xjj),final_matrix,'linear',0);%Ym contains a td-length vector for each point of x lateral distance for each element
                xjj=xjj+1;
            end
            %%%%%%%%%%% now calculating the weight (w) for each point:
            for adrs_z=1:length(Z)
                R=zeros(Lp,Lp);
                sum_Gp=zeros(Lp,td);
%                 SSS=zeros(P, Lp);
                for loopG = 1 : P
                    Gp=(shiftdim(Ym(adrs_z , : , loopG:loopG+Lp-1 ),1))'; % Converting from 3d to 2d.
                    R=R+ Gp*conj(Gp');
                    sum_Gp=sum_Gp+Gp; % (Gp.*repmat(WH', 1, td));
                    %*********** for subarray matrix-CF:
%                     SSS(loopG,:)= squeeze(Ym(adrs_z, tdd+1 , loopG:loopG+Lp-1));
                end
                R=R/P;
                R=R+(1/(20*Lp))*trace(R)*eye(Lp); %applying Diagonal Loading to R before calculating the weight.
                w=(R \ e) / (conj(e') * (R \ e));%  .* hann(Lp) ; % calculating the weight, w=(Lp x 1)
% %                 %%%%%%% CF %%%%%%%%%%%%%
                S=Ym(adrs_z, tdd+1 , :);
                CF=(abs(sum(S)))^2  / (N_elements* sum(abs(S))^2);
                %%%%%%% ESBMV:
%                 [w, Num]=func_ESBMV(R,w,0.5); % for Eigenspace-based MV: uncomment.

%                 %*********** for subarray matrix-CF:
%                 CF=(abs(sum(SSS))).^2  ./ (P* sum(abs(SSS)).^2) ;
%                 w=CF'; % .* hann(Lp);
%                 %*****************************
%  GCF=zeros(1,Lp);
% for ii=1:Lp
%     Spect=fft(SSS(:,ii));
%     GCF(ii)=sum(abs(Spect(1:M0+1 ,:)).^2) ./ sum(abs(Spect).^2);
% end
%                 [w, Num]=func_ESBMV(R,w,0.5); % for Eigenspace-based MV: uncomment.
%%%%%%%%%%%%%%%%%%%%%%
%                 S=sum(Ym(adrs_z, tdd+1 , :)); % output of DAS.
%                 out_power=abs(S).^2;
%                 Rp= sum((Ym(adrs_z, tdd+1 , :)-S).^2 )/N_elements;
%                 W_wiener=out_power* W_ESBMV /(out_power+ conj(W_ESBMV)'*Rp*W_ESBMV);
%                 w=W_wiener;
% % % % ************************************************************************
                B= CF*conj(w)' * (sum_Gp / P); % conj(CF) * (sum_Gp / P);  %   %conj(GCF) * (sum_Gp / P);   %CF * conj(w)' * (sum_Gp / P); 
                Beamformed_data(adrs_z,adrs_x)=B((td+1)/2); % taking the central point of the output vector.
            end
            adrs_x=adrs_x+1;
        end
    yy(1:size(Beamformed_data,1),1:size(Beamformed_data,2))=yy(1:size(Beamformed_data,1),1:size(Beamformed_data,2))+Beamformed_data;
    end
    
%%
% % % % Etot=sum(sum(Beamformed_data.^2));
% % % % Eout=sum(sum(Beamformed_data.^2 .* CM)); 
% % % % C= sqrt((Eout)/(Etot))
% % % xc=0/1000; % Place of cyst [m]
% % % zc=37/1000;
% % % RI=Relative_Intensity(Beamformed_data, X, Z, 2.5e-3, xc, zc)
    yy1=yy(1:size(Beamformed_data,1),1:size(Beamformed_data,2))/N_angles;  % averaging
%     env_dB=abs((Beamformed_data));
%     env_dB=20*log10(env_dB/max(max(env_dB)));
% % %     % ************************
% % %     [M N]=size(env);  % env_dB instead of J (J and env_dB have the same size)
% % %     sz=(M/15); % careful zainab!! the size is important
% % %     sx=(N/14);
% % %     Cst=env(round(5.24*sz) : round(8.76*sz), round(5.24*sx) : round(8.76*sx));
% % %     Bg=env(round(5.24*sz) : round(8.76*sz), 1 : round(3.52*sx));
% % %     CNR=(mean(Bg(:))-mean(Cst(:)))/sqrt(var(Bg(:))+var(Cst(:)))
    %% ************************
    env_dB=20*log10(abs(hilbert(yy1)));
%     env_dB=20*log10(abs(yy1));
    env_dB=env_dB-max(env_dB(:));
% env_dB=20*log10(abs(hilbert(Beamformed_data/max(Beamformed_data(:)))));
        %%%%%%%%%%%%%%%%%%% Display %%%%%%%%%%%%%%%%%%%%%%
    depth=z_start:L:z_stop;
    x=-image_width/2:dx:image_width/2;
 % %%%%%%%%%%%%%%%%%%%
    figure;
    imagesc(x*1000, depth*1000, env_dB, [-50 0])
    xlabel('Lateral distance [mm]')
    ylabel('Depth [mm]')
%     title(['CFMV: N=' num2str(N_elements) ' f0=' num2str(f0*1e-6) 'MHz,td=' num2str(td) ',Lp=' num2str(Lp)]);
    axis('image')
    colormap(gray(128));
    colorbar;
    toc
    return
        [MLW3dB FWHM FWHDR PSL DML] = calculate_performance_metrics(env_dB(1,:), dx, 50);
        Lat_Res=FWHM
        PSLL=PSL
        aaaa=env_dB(1,:);
        figure; plot(1:length(env_dB), env_dB(1,:)); 

    figure; 
%     imagesc(x*1000, depth*1000, env_dB .* CM , [-60 0])
    axis('image')
    colormap(gray(128));
    colorbar;
    %%
        [M N]=size(env_dB);  % env_dB instead of J (J and env_dB have the same size)
        Sz=(M/15); % careful zainab!! the size is important
        Sx=(N/14);

%         ff11=env_dB(round(0.5*Sz),:);
ff11=env_dB(2,:);
        figure;
        plot((1:N)/Sx-7, ff11); return
    % % %     plot((1:N)/Sx-5,env_dB(round(1*Sz),:) );   % for plotting lateral resolution at 30mm
    % % % % plot( (1:M)/Sz+25, env_dB(:,round(5*Sx)) );   % for plotting axial resoluton
    % % % grid on
    % % %     xlabel('Axial distance [mm]')
    % % %     ylabel('Intensity [dB]')
    %% %%%%%%%%%%%%%%%%%%%%%%% Lateral Resolution value:
%     [MLW3dB FWHM FWHDR PSL DML] = calculate_performance_metrics(env_dB(round(1*Sz),:), dx, 50);
%     Lat_Res(pitch_count,Lp_count)=FWHM;
% % %     Res_matrix(:,Lp_Counter, td_Counter)=env_dB(round(0.5*Sz),:);
%     end
% end
%% CF-weighted MV:
CFMV_Res=[-150.455540775757,-132.480193525813,-131.734165647260,-128.867002695880,-126.699879858749,-125.716357694119,-123.882323093151,-123.607957074482,-122.054737567774,-120.903387708200,-119.192531088113,-115.792471096543,-113.962185925897,-111.720199392767,-110.174142720077,-107.953706045789,-106.735797338078,-105.889318933235,-104.305614788871,-104.098307955889,-103.517684079224,-103.712538244676,-103.410026285474,-103.746772351108,-104.372099432333,-104.768314269582,-105.726337007976,-106.787958480557,-108.368996062858,-109.269638427103,-109.555089477872,-110.207703770537,-110.628395284590,-110.250857003551,-108.607573297870,-108.605401991342,-107.111357588652,-104.811192440817,-104.027631356087,-104.066713837184,-104.331580309813,-102.911042894677,-102.908066257802,-102.858613472715,-102.811342257189,-103.014973334167,-102.002867162282,-102.131530951476,-101.004637254403,-99.6088190722606,-98.9579308304983,-98.2773289696745,-96.7514176137678,-94.3167752436904,-93.6188064819725,-92.1941558842485,-91.3039012849634,-90.3768158780307,-89.5044582633951,-88.7517509955935,-88.3119530365485,-88.4262657290481,-88.3284492078531,-88.4046096243463,-88.7484863231490,-89.2760413555274,-89.5977795461270,-89.9466964567251,-90.6264449714354,-90.0709737443168,-89.3329968282209,-88.3964247972201,-88.4406490743221,-87.0164721623293,-86.6838875819312,-89.1409003318507,-93.9798866862257,-84.6368443094438,-77.5410221560626,-71.8231277160911,-67.1737471791679,-63.5681124701671,-60.5461640663937,-58.1887466772175,-56.1217789088472,-54.5083453234074,-53.3943627470779,-52.6882448838862,-52.3746922376095,-52.4985095117161,-53.1285303427421,-54.4275642276481,-56.5765751808968,-59.8742128706624,-65.6552285096930,-75.8699972703107,-74.0731884313681,-61.9296566116956,-52.6401924464999,-45.3708036151890,-39.5938596058437,-34.7323387229021,-30.4203903333948,-26.7697467457838,-23.3838043921747,-20.5825516262961,-18.1309342868918,-15.9009507130356,-14.1523542007436,-12.4245181634685,-10.6998874228054,-9.21452339580435,-7.98553513182458,-6.57806782759701,-5.00537257334804,-3.51348944020140,-2.29173924947247,-1.25057303222502,-0.498328885136289,-0.261596943279073,-0.0828624624458598,-0.0267882031681665,-0.0240548184972340,-0.0297449420404519,-0.165047037188515,-0.360296333244378,-0.809179690311908,-1.76107717473622,-2.80824964112185,-4.24301495024042,-5.77130720862078,-7.33099522675423,-8.57844374651995,-9.92300574307433,-11.5435488515831,-13.2732219305982,-14.9744921815065,-16.9420686276222,-19.3690621839086,-21.8849423247118,-25.0014218554853,-28.4782710535044,-32.4648616863736,-37.1003286116448,-42.2686490081671,-48.7603739049702,-56.8141901598957,-67.6939320573551,-78.9667392077952,-70.1712272566368,-62.3712131144501,-58.0656450423871,-55.4043466210006,-53.6846859631165,-52.7356527271544,-52.3860160477009,-52.4672549656668,-52.9959041998819,-53.8797658295180,-55.2551102654419,-57.0683367560804,-59.3057715795834,-61.9631947094758,-65.3171889267169,-69.3020525232228,-74.5674718911375,-80.5192636489875,-90.3509372152297,-91.0431896033750,-87.9055340414055,-87.0503655395173,-87.9320389641779,-87.9991083837948,-88.8445317196444,-89.5381697450487,-90.6764674277983,-90.7848174503276,-89.7982691598160,-89.3158724158924,-89.3649474260786,-88.7044448570987,-88.2044635208919,-88.3399118871141,-88.2563573162641,-88.4289557912337,-89.1671825309903,-89.8153884875036,-90.5834224570119,-91.6821177132106,-92.7406900799593,-94.0597013239640,-95.4904100325617,-97.4802860005311,-98.9850420893760,-99.6123947311723,-100.179766091609,-102.209409723112,-102.470581932453,-101.806309997774,-103.503502225962,-103.345338583691,-102.854969267790,-102.977860866459,-103.699609423851,-104.187437939823,-103.796256245648,-105.203164310762,-105.749924799690,-107.982579686420,-108.745510381180,-109.355757352006,-109.956500491009,-110.557317092139,-109.931261075943,-108.341957206869,-109.092888740238,-107.701951996494,-105.859670906198,-105.301111856242,-104.935722790188,-104.095942386338,-103.603515379939,-103.654593688606,-103.557477879005,-103.744677581484,-104.153627588976,-104.929161225063,-106.490282280113,-107.095893094392,-109.128957435336,-110.874673106587,-112.837035348909,-115.538382668982,-117.206601426975,-120.202853107652,-120.466519064415,-123.458102686170,-123.759186386561,-124.836930493993,-126.822491728693,-126.573729027223,-131.370310872095,-132.323881253431,-144.814496651909];
CF_weighted_Res=[-97.2302383173367,-96.2064970913247,-95.2386599799604,-94.6031932991407,-94.5133085888770,-94.0745210251420,-94.4608532211298,-95.0149784726446,-95.5061026316623,-96.9392166147083,-98.4412245741559,-100.588223568301,-102.712717412559,-104.542807358521,-105.172062784624,-104.226904809655,-102.790545367836,-101.680106265042,-100.398243636548,-99.4637507110368,-98.7482650618113,-98.2901121010829,-97.5441867008539,-97.1850529529164,-96.8461504202168,-96.2540906894952,-95.5248667767465,-94.9812690324600,-94.1115855063587,-93.7904268747271,-93.3956458123639,-92.3868805719955,-91.6896885376944,-91.1307289884632,-90.3290069609099,-89.5049641818304,-88.7313917953905,-88.1984144324926,-87.4171968362076,-86.9167523263581,-86.4354673400684,-85.8218777230520,-85.1812841104739,-84.5104240869012,-84.1724261089742,-83.7961426629725,-83.2693097435979,-83.0448208812125,-82.7311655976238,-82.5089075720242,-81.9758555601645,-81.6808189357859,-81.5185283065143,-81.0284740037536,-80.6341967461316,-79.8812770089417,-79.0229947028944,-77.6830996415212,-75.9263425825233,-74.1896499271692,-72.3110067131743,-70.1424057996479,-67.9161404904491,-65.7772611481126,-63.6774475313376,-61.6480644967942,-59.8241992218594,-58.2526837679943,-56.9157612900751,-55.6170228130031,-54.5292634708224,-53.5671870898394,-52.8539651971872,-52.3599661707389,-52.0496978698790,-52.0434306205679,-52.3344945719282,-52.8211992827422,-53.5513161943107,-54.7014014176049,-56.3533022469880,-58.4663256021822,-61.1842234694036,-65.0041215887734,-69.9517804398735,-74.2867927796411,-74.4870946683254,-73.1513883968561,-70.8068811825713,-66.5706341963453,-60.6770724022976,-54.0434590211539,-47.8041728257211,-42.2388461908913,-37.1808917583924,-32.5708883849829,-28.4825131852467,-24.8300117886312,-21.4890152366026,-18.4253139981951,-15.7054750391581,-13.2852554547526,-11.1576992313728,-9.29910510807605,-7.70600637418022,-6.37729609184470,-5.27033597100228,-4.32048611040716,-3.49999277058532,-2.81146166560688,-2.26554282123590,-1.83917830474400,-1.51198756427885,-1.24123586242507,-0.983802173873016,-0.737074539751347,-0.517678410667656,-0.343287747306249,-0.217015703928894,-0.130727885437352,-0.0673588577637361,-0.0345939647119735,-0.0264793296179846,-0.0471613594119731,-0.0943032511496540,-0.169526598674793,-0.272732806277247,-0.424452274987686,-0.620498385214944,-0.856177186157083,-1.11074148146560,-1.37239155840831,-1.66180585648607,-2.03475077281826,-2.52099367469788,-3.13206042590690,-3.88559265717851,-4.77466318005105,-5.79234037090470,-7.00432521847205,-8.45884454383685,-10.1743800406589,-12.1726088732112,-14.4480529277225,-16.9995037188068,-19.8990385754743,-23.1006779009152,-26.5850353759998,-30.4365395773646,-34.7796944280350,-39.6264191596056,-44.9032948069354,-50.8111174169877,-57.3477834742408,-63.7473054984988,-68.9294530596999,-72.1101487197287,-73.8132121289078,-74.8774924226601,-72.5101086327045,-67.4423620535940,-62.9545705746006,-59.7366430123338,-57.3791373192921,-55.4804642197419,-54.0636395969579,-53.1805735392787,-52.5484165776986,-52.1700707619318,-52.0020971682157,-52.1635910465976,-52.5687624475184,-53.1790100929995,-54.0021122534093,-55.0383035249151,-56.2441248449734,-57.5478933665334,-58.9955572961404,-60.7105901979728,-62.6328324148482,-64.6969039151623,-66.8260994932996,-68.9839853730301,-71.2745510683977,-73.2897361540999,-75.0444487508526,-76.7106418867925,-78.3470113100858,-79.5086809576243,-80.2810929316473,-80.8142556138575,-81.2070012504439,-81.5940400139003,-81.7636240391061,-82.4097971752097,-82.6493345749728,-82.8726282822682,-83.1712710855019,-83.3619721422393,-83.9733041500084,-84.2725063618422,-84.7825929328982,-85.4681959803185,-86.1693192186606,-86.6675136765080,-87.1903835378863,-88.0742825628389,-88.5508857053108,-88.9909068059371,-89.8922916142593,-90.6625573051052,-91.3646962584485,-92.0207163295802,-92.6638049660982,-93.3308182207134,-93.7948861307155,-94.5886406041984,-95.2758682305098,-95.9897821935040,-96.5874629713123,-96.9829043249162,-97.6677903748354,-98.1846653233833,-98.4690117084415,-99.0438992800279,-99.7693312418693,-100.761827879677,-102.306683148377,-103.466870580218,-104.644632118144,-105.474250529243,-103.520465769628,-101.911623388968,-99.2931062927949,-97.6459218759933,-96.4232142905328,-95.2341864794211,-94.6490168891899,-94.1301199665762,-94.2740122427917,-94.5411056706553,-94.8418669861284,-95.7139244895527,-96.3967350365313];
GCF_CF_weighted=[-102.336740823301,-101.466892452807,-100.614890957667,-100.038392519380,-99.9723614678723,-99.5486185542289,-99.9101861210017,-100.397215837879,-100.801543884572,-102.031594692956,-103.254906272017,-105.025490262056,-106.796477667866,-108.523063177636,-109.830104259817,-109.888049107453,-109.498818336170,-108.501069046572,-107.315762013351,-106.313406383797,-105.391915016832,-104.833490329700,-104.053441796226,-103.608111702897,-103.231829614700,-102.612191651425,-101.923471108425,-101.312513397202,-100.455057751446,-100.145965657332,-99.7565487212999,-98.7940974928162,-98.0875351998468,-97.5134921145126,-96.6719370534426,-95.8009498952326,-95.0144932661930,-94.4778318953130,-93.6991874510138,-93.1600011872732,-92.6656024932923,-92.0209887096927,-91.3419297922792,-90.6489510662600,-90.3061002929835,-89.9103783614132,-89.3631208878522,-89.1247463027282,-88.8015996843382,-88.5788252739209,-88.0194800148456,-87.7020624806951,-87.5290457090904,-87.0218138639289,-86.6232064810848,-85.8448024126909,-84.9688158084040,-83.6404888456069,-81.8914576068597,-80.1567170596267,-78.2791370716896,-76.1284401222159,-73.9035103876381,-71.7745900473030,-69.6796355845834,-67.6551356090386,-65.8344251336052,-64.2644480081400,-62.9307107015706,-61.6339613714446,-60.5469451959735,-59.5848110556250,-58.8666064683786,-58.3673172666876,-58.0495779900220,-58.0326337007562,-58.3115708197187,-58.7799196926746,-59.4918498990721,-60.6192704634784,-62.2464058375310,-64.3339637955493,-67.0278909826473,-70.8429001546731,-75.8692217608132,-80.6454092590923,-81.1570526024365,-79.7680708810451,-77.2265127740536,-72.6839741403089,-66.5877460439905,-59.8968959933401,-53.6774221465770,-48.1518973360430,-43.1358862567652,-38.5567275114192,-34.4785914107515,-30.8036225840470,-27.3950322060940,-24.2019432729500,-21.2816539551413,-18.5808262096843,-16.0798042555570,-13.7718813357135,-11.6573441960170,-9.76962823978556,-8.10036045231777,-6.61066838748241,-5.31110123227796,-4.20495330643655,-3.31024752834333,-2.60125505611950,-2.05102651032524,-1.61259024340058,-1.23530156682568,-0.901128742929473,-0.619037309565044,-0.400353318112821,-0.240048975517766,-0.132784638507815,-0.0625609692727949,-0.0333919411439183,-0.0271565783610299,-0.0434127411709824,-0.0914915178371984,-0.180246003288460,-0.311432373701564,-0.501623210554158,-0.749708651647268,-1.05993674457579,-1.41644116522929,-1.82024771643972,-2.30367575364392,-2.92895870651898,-3.72909522537708,-4.72386217871724,-5.92532001437860,-7.32524893082314,-8.90085255402244,-10.6739500791996,-12.6743662576990,-14.8767600121464,-17.2898454633983,-19.8942533543079,-22.6840789058648,-25.7487156316840,-29.0473865791007,-32.5758353821515,-36.4324482792438,-40.7540695183169,-45.5625282042456,-50.7964865221491,-56.6701945083320,-63.2127458326810,-69.7350142446192,-75.1792667990427,-78.6210873607060,-80.4368947438966,-81.4623181228157,-78.5829691000779,-73.2943435872174,-68.7879752394783,-65.5883642873904,-63.2580939278953,-61.3837468777498,-59.9912149524870,-59.1282903562168,-58.5137289741007,-58.1509385574484,-57.9933747042298,-58.1643404598467,-58.5747008427724,-59.1907121986022,-60.0163386405051,-61.0541726547758,-62.2613951963447,-63.5640720815709,-65.0113130863286,-66.7255180090492,-68.6430804486314,-70.7035495002278,-72.8263698359876,-74.9815253678321,-77.2677896842600,-79.2694401079049,-81.0259310728299,-82.6927947759535,-84.3266791596010,-85.4727085900252,-86.2711167775080,-86.8048116107005,-87.1991823498128,-87.6022707699755,-87.8011607556185,-88.4636866310101,-88.7196609872726,-88.9439910478479,-89.2673757834075,-89.4574059405890,-90.1041592130948,-90.4203577785979,-90.9498462270993,-91.6645670200522,-92.3984559944740,-92.8989753850484,-93.4533292120304,-94.3798446392268,-94.8644836114993,-95.3043445964877,-96.2464321691651,-97.0502918169199,-97.7679433330926,-98.4447174893612,-99.1098310527003,-99.7472420261685,-100.209078659881,-100.930736890778,-101.686625964999,-102.406268945338,-102.991982310152,-103.453715121339,-104.066710717682,-104.637064108242,-105.076967433863,-105.828326764748,-106.760469947016,-107.815305586354,-109.111417250186,-109.763945132419,-110.252906675959,-109.691186951742,-107.387359983638,-106.043476380400,-103.941155107601,-102.597273713329,-101.589462449069,-100.553639809113,-100.050277671818,-99.5716142901540,-99.7135849847097,-99.9600083022780,-100.213329226458,-100.996775272643,-101.570815901216];
figure ;
p=(1:N)/Sx-2.5;
plot (p,CFMV_Res, p, CF_weighted_Res, p , GCF_CF_weighted);
grid on
title ('CF-weighted ESBMV');
xlabel('Lateral Distance [mm]');
ylabel('Amplitude [dB]')
legend('CFMV','CF-weighted MV', 'CF-weighted MV+GCF');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
CFWMV_latR=[-27.3811085954037,-30.1763353838987,-28.3514545473623,-24.2082429342240,-21.7790228801886,-20.6409714418777,-20.6438036753665,-21.9122976166075,-24.6128040349860,-28.7686369799782,-34.4811739814608,-37.6886063054146,-30.1811171197214,-25.7029965682053,-24.1873801252417,-25.5439287496055,-31.2856008391280,-38.0022009959726,-25.1872118869155,-20.6373354252353,-17.8596801311762,-15.0634626248942,-11.8971529451352,-9.50976307531442,-8.25579080870011,-7.72990583682514,-7.97209355719258,-8.90961708525629,-10.3723613389014,-12.4058314070396,-15.0189454883827,-17.8390340384288,-20.2521505168576,-22.1758990511487,-23.9859139298461,-25.8447504205525,-27.8703090599459,-29.9220039201012,-31.9247289196753,-33.8422565873710,-35.6096244931934,-37.3830706534193,-39.7201483096678,-42.6622227837835,-46.6761173770802,-52.1876150422998,-59.1063940960369,-60.2547114379772,-56.6479063218600,-53.1570826605097,-49.7972887742127,-46.5360564040251,-43.8611138817241,-41.7206631333369,-40.6350171559254,-40.1756517802469,-40.5172452530521,-41.5463023223464,-43.0528206462604,-44.6429377690058,-46.7346714323380,-49.0419013838529,-50.3942698375657,-43.0003137466069,-36.5843234343818,-32.4486400928424,-29.6863353846402,-27.9630789482801,-27.1275325241564,-27.0269221247893,-27.8510746836014,-29.4147986110740,-31.7378392118266,-34.9354668368387,-39.3675180255083,-45.8002821257806,-57.4989700777699,-60.8417091353239,-56.1548087219068,-57.1206623116137,-61.2772474499885,-66.7521615946423,-68.9905785241378,-68.8158251980531,-69.0074153242335,-72.5668533237050,-70.6419966654220,-62.9298667998840,-58.3236970182165,-56.1476923447432,-54.9491442430474,-54.7585991712390,-54.9578984206520,-55.4520740687474,-56.7780471509831,-58.7315200840733,-61.5635851769604,-65.2935527405154,-70.0693924905109,-75.1627361340420,-80.7252175753399,-93.0032735474899,-83.9287204243848,-77.0723699554804,-73.1764532104167,-70.7324552902859,-68.9882678989148,-68.6253126802864,-69.2189401656236,-70.2804761159691,-73.4152049669220,-78.9259248466239,-88.7772680136694,-80.0032865858958,-77.4191346661629,-76.9785428562312,-78.5490196519035,-80.1961897719405,-79.5012589448547,-76.2499734435343,-72.3750234892670,-69.5646914189723,-67.8004050150353,-67.3080361847022,-67.5741532880620,-69.5705976264893,-73.5839289729841,-80.2152456424097,-74.8574779336004,-70.0397451100465,-67.9132377340263,-66.6584645822681,-65.8769491559808,-64.1304235358639,-60.9365361845675,-57.6141173589888,-54.3000888506047,-51.7726554827731,-50.4154524789842,-49.8967437665437,-50.6817054514809,-52.6990762072205,-56.7592794631109,-60.8099309762092,-56.4910535302403,-52.6864497524145,-50.4668477209045,-48.9102175997282,-47.6258592974570,-46.5717841482312,-45.9118416469541,-45.7206881741838,-44.6955606059173,-40.4597273371771,-34.9895772033973,-30.1797734480312,-26.2636215825769,-22.9143490372708,-20.0345062625229,-17.6511884894996,-15.6681701060510,-14.1123783786287,-12.9974135290380,-12.2020522959685,-11.7382315286970,-11.5416392611610,-11.5477998953749,-11.8558178591715,-12.3915261060594,-13.5393757858526,-15.4658455503795,-17.9517197892024,-21.0461479481313,-24.2998280821492,-27.7354234679114,-32.0345262247550,-41.8181883985251,-40.7625491690225,-32.2515346781315,-29.9260306851562,-30.3222784677184,-32.1909573218240,-34.2047900059444,-35.6332596401834,-36.4499783995345,-37.6662448592830,-38.7437238762096,-38.8900287550292,-37.6271445126450,-35.7363108375035,-33.5128957168908,-31.3733075468911,-29.5634592521140,-28.2238960176943,-27.4008204049370,-27.3281033237442,-27.8902255486935,-29.2818341356369,-31.2211101195875,-32.8622045691234,-32.0282670228199,-30.3053093571077,-29.1879243049193,-29.5112408692268,-31.0872210819328,-30.1420893720380,-25.0796370823905,-20.8560072728754,-17.9066078491738,-16.4993269044073,-17.2306472320544,-20.5950577203163,-25.5372866130585,-25.8503225363182,-23.0762270053933,-20.2997474469233,-17.8058951371841,-16.1001240584837,-15.7146947793410,-17.0683989233251,-20.3835328365659,-26.2131509696920,-30.2129470214083,-29.5707460038350,-26.6836131160850,-22.8011209100843,-19.3187143638864,-17.1523894219883];
CFMV_latR=[-23.0511132784384,-38.6802790795330,-28.1017308887725,-27.2629959591620,-24.6753043793865,-22.0710718732046,-20.8219164733816,-21.1630305348355,-23.2102915854779,-27.0334844971136,-32.8402882984775,-43.6029914697877,-41.5795961170192,-37.5394558977615,-36.7329500637382,-38.4447607567348,-43.7027272824408,-61.5855514854102,-51.6368609769564,-44.4170529139189,-35.4736434930485,-28.0582293034362,-22.7053264465471,-19.2495005126132,-17.2435593904543,-16.3892616750792,-16.5266359368721,-17.6897498397192,-19.4440393589684,-21.3021905929906,-22.7661459667539,-23.9769393862015,-25.3825160574849,-27.1549261350103,-29.4516917030800,-32.5116829863895,-36.3011717938112,-40.4279622255954,-45.0106976905402,-49.6367287642178,-54.7545008874035,-59.5421955203260,-62.9836119437542,-63.5316774237829,-62.3811689243040,-60.9521785754409,-59.8529185148782,-58.8408570073462,-57.0752803708697,-54.3863926495851,-51.6149650691957,-48.8169364024106,-46.4449262071929,-44.5743693856693,-43.9842219493769,-44.8417684621690,-47.3603723483508,-52.0699020277054,-59.9808685764932,-74.1925214780966,-58.4206191506874,-53.6208201669801,-51.9596779916693,-43.7798169876623,-36.7921347835384,-32.5653280418988,-30.0314692482773,-28.6621196822070,-28.3394646679405,-29.0073278169227,-30.8369886753316,-33.6002728362725,-37.2120791115732,-41.5407282872536,-46.8976314547563,-52.3301323987013,-57.6178496353576,-63.2840914899951,-70.1039245194831,-81.8678361658447,-86.2399443680305,-79.4263711923585,-75.9778538576618,-73.5296225236139,-72.6534244170976,-74.1848387889913,-77.8497598902196,-83.8972088366521,-91.7141756756725,-93.2297608437778,-90.5672019224630,-90.2223067923727,-95.7521977507378,-103.194755022068,-91.8724131592900,-90.9831589853583,-94.1747631734779,-102.651430658939,-107.362270439415,-104.483697648492,-100.649843706736,-95.4894725497085,-94.2187881841736,-94.5770607683122,-97.9362893785323,-104.812694251410,-112.164104460069,-104.327159142222,-102.898898644484,-103.517603352668,-104.091914584498,-101.188115063213,-99.1154237604878,-100.244376374285,-100.757645000405,-100.750495136396,-100.402854206806,-104.267996402838,-108.775604809372,-98.2131722251300,-94.0310874692518,-89.7494533849011,-88.1703883273812,-88.0084230376688,-89.3335273710580,-93.5951044443094,-99.6262930423877,-98.3477510101041,-95.5666226063348,-89.0034703261562,-82.4363866390750,-76.4380062914261,-72.3054179704629,-70.1409298797635,-70.0116251252301,-72.3814363275109,-78.7586960228593,-96.6203760611573,-81.7720202448323,-70.7246811779682,-65.0496781262984,-62.6779956613006,-63.6745755661414,-68.9761341600029,-76.2915095881126,-68.1916440398541,-68.0830690540171,-74.8770997365071,-93.1797938431571,-70.0295748574030,-58.5813681395228,-52.7926416304650,-50.6077516848546,-47.9348929758264,-42.9141362176877,-37.8578941504149,-32.9792838217127,-28.0118392901713,-23.2642628077013,-19.1226094529608,-15.8030399678069,-13.5021713248892,-12.0864200505540,-11.3284996718424,-11.2828019538043,-11.7959038935786,-12.6219618192175,-13.7394024597079,-15.0994481958563,-16.8471892116378,-19.2225100835545,-22.6180055231109,-27.2415373285110,-33.9878963129825,-44.6806285127013,-58.4608131858432,-63.1742932309125,-61.3299326613591,-47.0632512641585,-38.5564066231175,-33.1867694806548,-29.5562584410711,-27.2739831986316,-25.8811416113959,-25.3351868187995,-25.5177853898006,-25.9487276218866,-26.4322802360369,-26.3961663366370,-25.9751499218265,-25.2170518865201,-24.8008943828619,-24.9262279577459,-25.9302267110765,-28.0313381702463,-31.5584327193239,-35.5613834330979,-37.0076738511828,-35.8014194222085,-33.2365563717216,-29.5371039893504,-25.0470560897130,-19.8930439097175,-14.3704217342022,-9.62237052004747,-6.84573546300834,-5.30292903252001,-4.41255196508752,-4.17479212500189,-4.82706674996877,-6.78297062048870,-9.73856974315839,-12.4961358378667,-17.9532153831192,-30.2941904048005,-13.5468099933251,-8.77831306954278,-7.08295175570407,-7.90551531100203,-11.7130421962228,-18.2463620537697,-28.4759343524223,-48.4471142366709,-47.4257539819554,-28.9848348881231,-19.5288293836628,-14.8557789229516,-12.7505567627945];
p=(1:N)/Sx-7;
figure;
plot (p,CFWMV_latR, p, CFMV_latR);
grid on
% try2:
CFWMV=[-30.8800510758387,-39.3117837004720,-22.9074718762316,-15.5848693666250,-12.2012393282022,-10.7493235974628,-10.6482980965721,-11.7699910007412,-13.9430494129045,-17.6699251776124,-21.9134304999498,-26.6239330209626,-30.2199955243119,-32.5597146339163,-34.0136194756846,-34.9257800601978,-34.8051894171678,-33.6438362638010,-32.1991513012435,-30.8433791552889,-28.5720378246045,-25.3230058807854,-21.8228660944458,-18.9118054152516,-16.7398294432014,-15.2132140479002,-13.9546393297591,-13.1141669301531,-12.6169986504133,-12.4825888491538,-12.8001138894958,-13.1997118513375,-14.0363858517854,-15.3786103984270,-17.1012495454242,-19.2346866618527,-21.9104517554633,-24.6243834194474,-27.5774107549593,-30.9646550728816,-35.0141458491512,-39.7131538265346,-45.1969793915957,-50.3937624492795,-53.5080970367443,-54.3715805674761,-53.0458411662374,-50.6802609472131,-48.1807287846966,-45.3934795790631,-42.3462303248236,-39.5160503493626,-37.1248252408191,-35.4358457700041,-34.7195417793093,-35.0522654613293,-37.0184368015225,-40.5625390648778,-46.3784273768424,-53.4515903238006,-60.7079458485084,-54.9798047388313,-40.3226660345556,-30.5127392068954,-23.9359657797123,-19.3863220021817,-16.2671747260795,-14.1782589999315,-12.9742438102963,-12.4491071260321,-12.8083616531900,-13.9657813437030,-15.8054349940734,-18.4901189071140,-22.0465371607705,-26.2675383727363,-31.5394792231778,-38.0431286945333,-45.7231887943370,-54.5569603546614,-61.1991190571146,-64.0070318182989,-63.2223769588614,-58.7291184097293,-55.2125801665898,-54.0311069029369,-54.5373853416216,-56.4691278658659,-59.4948922823415,-61.4114755441826,-61.1420039777014,-60.2086409709458,-60.7462227841397,-61.9422970849353,-63.2993058760618,-64.0497501874186,-63.6664990197822,-64.0996645852857,-65.3472121366848,-67.3898441694541,-71.2580359677240,-82.0084336463173,-77.3627995215014,-70.0963327594985,-66.2424029469694,-64.8189184878536,-63.9677430001811,-63.6652092776296,-64.6183758043674,-65.8712359167861,-69.5023160260234,-78.5857889988495,-78.0077613562340,-71.1114247037906,-69.9363638586734,-70.3522927385856,-71.6813859249089,-72.0157921466855,-68.9507545913442,-65.3938486078566,-62.1429315427467,-59.9700619651949,-58.5830648011212,-58.1971895778335,-58.2553449274842,-59.4384407708182,-60.7968958749523,-60.8511872824033,-59.6867744972203,-58.3593049131476,-58.1020742809144,-58.4117870204407,-59.2219340992230,-58.5996168513441,-53.6432194315457,-48.5358887953960,-44.2489895629947,-41.1175756730327,-39.6431128550944,-39.2606502022588,-40.0668813456243,-42.1398770561470,-46.2860604222672,-51.1779656343565,-48.8930118663100,-46.2190723291707,-45.3020290374831,-45.3079846395551,-45.3744351739949,-43.6007748516557,-41.1289029580141,-39.0109003005417,-37.7555789158358,-34.1351567096808,-27.6891289616526,-21.6261555420584,-16.5486144104216,-12.4292877468953,-9.01600906106637,-6.26744142913697,-4.12750406667271,-2.77946440197246,-1.95275857970114,-1.65326306529789,-1.75891397294720,-2.24352528068914,-3.13509678725723,-4.61404812383887,-6.44323414815221,-8.96372446888802,-12.0572445800577,-15.7116683003185,-20.1995208374099,-25.4727786754195,-31.9732439963619,-40.1829794711940,-50.9904430093557,-49.7895634270797,-42.4408633932893,-37.8719820751732,-35.3255704859656,-34.6452046974177,-34.6218897028631,-33.9631368093620,-32.8468495310317,-31.1081800210866,-28.7744175280910,-26.2813249449281,-23.8946168379169,-21.5140667595777,-19.3382725110582,-17.4033004444265,-15.8962214190619,-14.9376379002093,-14.5949538222513,-14.9391488377806,-15.9741609592567,-18.2362051457799,-21.9382147150837,-27.9325302355824,-38.4233188096601,-52.5396577460066,-47.9382523598704,-36.5541346505693,-26.3527275532986,-19.3986360664196,-14.3959675614392,-10.7053213408322,-8.16617565210589,-7.09398155455551,-8.27207806538928,-11.9776577003053,-17.2683749626369,-17.4885092511186,-11.2504646462728,-6.54060280391110,-4.11938908027139,-3.27227089082743,-3.78162804330941,-5.79675159644830,-10.3183783242692,-19.3990065824012,-33.6891711149144,-21.2484225852243,-12.4432590990226,-8.10328111546232,-5.64680915166588,-4.57622838851319];
CFMV=[-27.5788916951461,-49.6071754551191,-40.8268230453897,-24.6630627050884,-18.6141532951116,-15.5184297488686,-14.0646731718589,-13.9774870780990,-15.0746176983746,-17.8243615599409,-21.6262499021162,-27.1659364030918,-34.7634167680470,-45.7120414611409,-65.4191674375636,-67.2998756852662,-67.2232121088449,-53.4883510378234,-51.4715405213727,-40.3915222692957,-29.7789809398216,-23.0286152946570,-18.5145120950336,-15.5250447281929,-13.6498894760309,-12.5765340821972,-12.0380018660074,-12.0166891642513,-12.3712710188129,-13.0763811926973,-14.1153563313685,-15.2740212589151,-16.5850846536645,-18.3382265024133,-20.4024450285987,-23.0308955855609,-26.5297438448309,-30.4704968182946,-35.0832846243738,-40.1948657477847,-45.4837099615185,-50.0164344281415,-53.7230200904767,-55.6137812467337,-56.1255377807992,-56.0262463312948,-54.9837981654365,-53.3031847686432,-51.3350390891363,-48.8214579802668,-46.1794313830618,-43.4004518855400,-40.7996531087861,-38.6237988507104,-37.4974285251192,-37.6832482237103,-39.8181271904715,-43.9840185208581,-51.3309436685371,-63.6345187088957,-77.2718953072977,-60.7848938458319,-42.9946894550445,-32.4354570448322,-25.6926714923301,-21.2898648757378,-18.5074002493318,-16.7511588168132,-15.9176324786499,-15.7757962286267,-16.5584194256260,-17.9922519882034,-20.0052305402616,-22.8966000221692,-26.7759213603847,-31.4406131852933,-37.3257811847674,-44.7051047993923,-54.1514704392297,-68.3414201390997,-91.9365621341335,-98.8037915738721,-78.0164476874303,-68.7326718163004,-64.6606514316384,-64.5369327689683,-66.8357159010691,-71.5121597346414,-79.5839458760402,-91.2309193713875,-107.063719608197,-123.716338030939,-96.2078778644490,-87.5831019616929,-84.2315717687931,-85.3646849242398,-89.9133497415729,-102.517450369933,-102.460543166506,-98.3397974508916,-95.1085797046460,-93.1581128504522,-93.8919006853739,-96.0540272159034,-96.5371014380359,-98.6937737425769,-102.785152283018,-112.259118378352,-114.223025946356,-103.725854818837,-95.4114185952649,-91.7762558214918,-90.8883070541554,-93.6730231355879,-101.211031279330,-103.929622990108,-97.2113222169704,-98.7064983548719,-123.978749443969,-93.1173810429356,-84.9428786628239,-81.4256852187853,-80.2501241654705,-80.5457981059642,-82.4933090445574,-88.5696920167358,-100.460398138499,-93.1022559426381,-102.088169216022,-81.5192040410142,-73.2285901975940,-68.1150758204144,-65.0655441196032,-63.6235887500213,-63.8154595871787,-66.3874054852931,-74.2630754652141,-78.3558277022314,-66.6173300606764,-60.2561845966059,-56.4794254517035,-55.2103762089754,-56.8214629265005,-61.5125997459018,-64.5444000212671,-61.1769159082896,-62.1621898801892,-72.3876394417686,-73.6514484004047,-60.9930201060572,-51.6020345046669,-47.4672369052472,-50.4632185913520,-54.1096468226833,-42.8136693981529,-47.8300514472061,-33.4050311009360,-22.2997367598422,-15.6012013600252,-10.8704248066581,-7.51792156051943,-5.59313182479286,-4.45365966934298,-4.08576101086055,-4.35538785388917,-4.97973690112872,-5.76165944343256,-6.87295409665683,-8.09851176754421,-9.72337491610290,-12.0000649243040,-15.0051970814452,-19.4478133168085,-25.8903205421151,-35.9673093767299,-54.3090398882921,-95.3605413858529,-55.1939080562522,-40.7705982773034,-34.0482238589934,-30.8480948989070,-30.1329486773493,-30.9177150635253,-32.0196389074661,-32.8588794418032,-32.3465367408677,-30.0306075357146,-26.7719247660434,-23.5019938474234,-20.5988260319695,-18.5500524482283,-17.3156161459453,-17.0143592668801,-17.8667140405577,-19.7360915654933,-22.6946809167574,-26.4224215790853,-31.0967348681152,-36.1936565679907,-41.7951623116040,-48.0895798624044,-51.4030269267453,-51.3418230890187,-57.6607739660587,-63.4279907593163,-42.9700124612773,-28.8782493466177,-20.2579239363045,-14.7928987573664,-12.1343100597552,-12.5498365042872,-16.6837151221557,-25.5303489265855,-26.5595595794491,-14.6294509281914,-8.13339475934697,-5.20564263871347,-4.37178888458300,-5.13006481394427,-7.62257889073976,-12.9278464851816,-23.8576910591074,-54.1871373527013,-26.0400585958830,-15.6770980403288,-10.1191405125147,-6.76193279935046,-5.45678021426079];
p=(1:N)/Sx-7;
figure;
plot (p,CFWMV, p, CFMV);
grid on