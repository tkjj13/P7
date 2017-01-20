%==========================================================================
%
%  By Gennady Zilberman, Ben-Gurion University of the Negev
%  genazilberman@yahoo.com
%
%   This software performs Forward Channel simulation of IS-95 standrd
%   system. 
%
%   The simulation environment: Static AWGN channel.
%   Link Rate = 9600 KBps
%
%==========================================================================


global Zi Zq Zs show R Gi Gq

clear j;
show = 0;
SD = 0;         % ---- Soft/Hard Decision Receiver

%==========================================================================
% ---------------------------------------- MAIN SIMULATION DEFINITIONS --------------------------------
%==========================================================================
BitRate = 9600;
ChipRate = 1228800; 
N = 184;  		         % 9.6 KBps rate -> 184 netto data bits in each 20 msec packet
MFType = 1;		     % Matched filter type - Raised Cosine
R = 5;			            % Analog Signal Simulation rate

% ------------------------ Viterbi Polynom -------------------
G_Vit = [1 1 1 1 0 1 0 1 1; 1 0 1 1 1 0 0 0 1];
K = size(G_Vit, 2); 		% number of cheap per data bit
L = size(G_Vit, 1); 		% number of cheap per data bit
  
% ------------------------ Walsh Matrix -------------------
WLen = 64;
Walsh = reshape([1;0]*ones(1, WLen/2), WLen , 1);
%Walsh = zeros(WLen ,1);

%------- Spreading PN polynomials ---------
%Gi = [ 1 0 1 0 0 0 1 1 1 0 1 0 0 0 0 1]';
%Gq = [ 1 0 0 1 1 1 0 0 0 1 1 1 1 0 0 1]';

Gi_ind = [15, 13, 9, 8, 7, 5, 0]';
Gq_ind = [15, 12, 11, 10, 6, 5, 4, 3, 0]';

Gi = zeros(16, 1);
Gi(16-Gi_ind) = ones(size(Gi_ind));
Zi = [zeros(length(Gi)-1, 1); 1];     % Initial State for I-Channel PN generator

Gq = zeros(16, 1);
Gq(16-Gq_ind) = ones(size(Gq_ind));
Zq = [zeros(length(Gq)-1, 1); 1]; 	  % Initial State for Q-Channel PN generator

%------- Scrambler Generation Polynom ---------
Gs_ind = [42, 35, 33, 31, 27, 26, 25, 22, 21, 19, 18, 17, 16, 10, 7, 6, 5, 3, 2, 1, 0]';
Gs = zeros(43, 1);
Gs(43-Gs_ind) = ones(size(Gs_ind));
Zs = [zeros(length(Gs)-1, 1); 1];  		% Initial State for Long Sequence Generator


%------- AWGN definitions -------------
EbEc = 10*log10(ChipRate/BitRate); 		% Bit/Chip rate loss
EbEcVit = 10*log10(L);
EbNo = [-2 : 0.5 : 3]; 		% measured EbNo range (dB) 

%==========================================================================
% ----------------------------------------------- MAIN SIMULATION LOOP ---------------------------------
%==========================================================================
ErrorsB = []; ErrorsC = []; NN = [];
if (SD == 1)
   fprintf('\n SOFT Decision Viterbi Decoder\n\n');
else
   fprintf('\n HARD Decision Viterbi Decoder\n\n');
end

% JHM DEFINITIONS/CHANGES 
% 
ErrB_N = 30 % ORIGINAL VALUE 300
iter_N = 150 % ORIGINAL VALUE 150

% WCS7 Amplifier data
data = dlmread('../../DATA01.R1',',',9,0);
S11 = data(1:end/4,:);
S21 = data(end/4+1:end/2,:);
S12 = data(end/2+1:3*end/4,:);
S22 = data(3*end/4+1:end,:);
p = linspace(-15,10,201);


for i=1:length(EbNo)
   fprintf('\nProcessing %1.1f (dB)', EbNo(i));
   iter = 0;	ErrB = 0; ErrC = 0;
   while ( ErrB < ErrB_N ) & ( iter < iter_N ) % JHM INTRODUCED VARIABLES
      drawnow;
      
      %----------------------------------- Data Transmit ---------------------------------
      TxData = (randn(N, 1)>0);          %- Gamble the Tx data
      
      [TxChips, Scrambler] = PacketBuilder(TxData, G_Vit, Gs);		% 19.2 Kcps rate
      [x PN MF] = Modulator(TxChips, MFType, Walsh);  % 1.2288 Mcps rate
      
      %-------------------------------- PA model ------------------------------------
      % JHM: 'x' is the transmitter signal and your PA model should be applied to this signal
      % NOTE: REMEMBER TO SCALE THE POWER LEVEL OF THE SIGNAL BEFORE ADDING
      % YOUR MODEL!
      % S_TX = x;

    % transform from cart 2 pol
        [PolPhase PolAmp] = cart2pol(real(x),imag(x));

    % calculate power of signal
        Power = mag2db((PolAmp.^2)/2);

    % find amp stats for power

        for k = 1:length(Power)
            index = 1;
            while Power(k) > p(index)
                index = index+1;
            end
            alpha(k) = abs(S21(index,1)+1i*S21(index,2));
            phase(k) = angle(S21(index,1)+1i*S21(index,2));
        end


    % calculate end power and phase
        AmpO = PolAmp.*transpose(alpha);
        PhaseO = PolPhase+transpose(phase) - repmat(angle(S21(1,1)+1i*S21(1,2)),size(transpose(phase)));;
    % trasform from pol 2 cart

    [sOI sOQ] = pol2cart(PhaseO,AmpO);
      S_TX = sOI+1i*sOQ;
      %-------------------------------- AWGN Channel ------------------------------------
      % JHM: Here is where you add noise. I have left the code from the original file here for your
      % reference. As you can see the noise is not added correctly. 
      % NOTE: REMEMBER TO SCALE SIGNAL AND NOISE CORRECTLY. 
      %
      
      noise = 1/sqrt(2)*( randn(size(S_TX)) + j*randn(size(S_TX)))*sqrt(R/2)*10^(-(EbNo(i) - EbEc)/20);
      noiseAmp = randn(size(S_TX));
      noisePhase = randn(size(S_TX));
      Power = mean(abs(S21(:,1)+1i*S21(:,2)));
      noise = Power/sqrt(2)*sqrt(R/2)*(noiseAmp.*cos(noisePhase)+1i*noiseAmp.*sin(noisePhase))*...
          10^(-(EbNo(i) - EbEc)/20);
      
      %
      r = S_TX+noise;
      
      %------------------------------------ Receiver ---------------------------------
      RxSD = Demodulator(r, PN, MF, Walsh); 		% Soft Decision   19.2 Kcps
      RxHD = (RxSD>0); 		% Define Hard Decisions received chips
      if (SD)
        % RxSD = sign(RxSD); 		% SOFT Decision = HARD Decision
        % RxSD = sign(RxSD).*(log2(abs(RxSD)+1)); 		% SOFT decision = 3 bits representation
        [RxData Metric]= ReceiverSD(RxSD, G_Vit, Scrambler); 		% Get Received Data 9.6 KBps (Soft Decision)
      else
        [RxData Metric]= ReceiverHD(RxHD, G_Vit, Scrambler); 		% Get Received Data 9.6 KBps (Hard Decision)
      end
      
      
      if(show)
         subplot(311); plot(RxSD, '-o'); title('Soft Decisions');
         subplot(312); plot(xor(TxChips, RxHD), '-o'); title('Chip Errors');
         subplot(313); plot(xor(TxData, RxData), '-o'); 
         title(['Data Bit Errors. Metric = ', num2str(Metric)]);
         pause;
      end        
       
      if(mod(iter, 50)==0)
         fprintf('.');
         save TempResults ErrB ErrC N iter
      end
      
      ErrB = ErrB + sum(xor(RxData, TxData));  % Data bits Error
      ErrC = ErrC + sum(xor(RxHD, TxChips));  	  % Chip Error (before Viterbi)
      iter = iter+ 1;
   end
   %---------------------------------------- Save the data of current  iteration ---------------------
   ErrorsB = [ErrorsB; ErrB];
   ErrorsC = [ErrorsC; ErrC];
   NN = [NN; N*iter];
   save SimData *
end

%---------------------------------------- Calculate Error's probability ----------------------------
PerrB = ErrorsB./NN;
%PerrB1 = ErrorsB1./NN1;
PerrC = ErrorsC./NN;
Pbpsk= 1/2*erfc(sqrt(10.^(EbNo/10)));
PcVit= 1/2*erfc(sqrt(10.^((EbNo-EbEcVit)/10)));
Pc =   1/2*erfc(sqrt(10.^((EbNo-EbEc)/10)));

%% ---------------------------------------- Show simulation Results -------- ---------------------
figure; 
semilogy(EbNo(1:length(PerrB)), PerrB, 'r-o'); hold on;
%semilogy(EbNo(1:length(PerrB1)), PerrB1, 'k-o'); hold on;
semilogy(EbNo(1:length(PerrC)), PerrC, 'm-o'); grid on;
semilogy(EbNo, Pbpsk, 'b-.x'); xlabel('EbNo (dB)');   
%semilogy(EbNo, PcVit, 'k-.x'); ylabel('BER');
semilogy(EbNo, Pc, 'g-.x'); 

title('Error probability as function of EbNo');

legend('Pb of System (HD)', 'Pb of System (SD)', 'Pc before Viterbi of System', ...
   'Pb of BPSK with no Viterbi (theory)', 'Pc on Receiver (theory)');



legend('Pb of System', 'Pc before Viterbi of System', ...
   'Pb of BPSK with no Viterbi (theory)', 'Pc before Viterbi (theory)', 'Pc on Receiver (theory)');
