close all;clear;clc ; 
%###############################
Ts = 1 ;
step=0.01;
t = 0:step:(Ts-step) ;
t_length = length(t);
Rs = 1./Ts ;
number_of_bits = 10000 ;
%##############################
sent_signal=zeros(1,t_length*number_of_bits);
received_signal=zeros(1,t_length*number_of_bits);
ss=1;
%###############################
%BFSK Signals
f1 = 1 ;
f2 = 2 ;
rand_ph=randn(1,1)*pi;

signal_1 = cos(2.*pi.*f1.*t+rand_ph) ;
signal_11 = signal_1';  
signal_2 = cos(2.*pi.*f2.*t+rand_ph) ;
signal_22 = signal_2';   
%###############################
cos_1=cos(2.*pi.*f1.*t) ;
cos_11=cos_1';
sin_1=sin(2.*pi.*f1.*t) ;
sin_11=sin_1';
cos_2=cos(2.*pi.*f2.*t) ;
cos_22=cos_2';
sin_2=sin(2.*pi.*f2.*t) ;
sin_22=sin_2';
%###############################
k = 0 ;
Eb_N0_dB = 0:12;
BER_simulation = zeros(1,length(Eb_N0_dB)) ;
BER_theorical = zeros(1,length(Eb_N0_dB)) ;
BER_theorical_coherent = zeros(1,length(Eb_N0_dB)) ;
%###############################
%TX
tx_bit_stream = randi([0 1],1,number_of_bits);
tx_signal = zeros(t_length,number_of_bits) ;
    
    for ii = 1:1:number_of_bits
        for m = 1:t_length
            if(tx_bit_stream(1,ii) == 1)
                tx_signal(m,ii) = signal_22(m,1) ;
            elseif(tx_bit_stream(1,ii) == 0)
                tx_signal(m,ii) = signal_11(m,1) ;
            end
        end
    end
%###############################

for Eb_N0_dB = 0:12
    Eb=sum(signal_2.^2);
    Eb_N0 = 10.^(Eb_N0_dB./10) ;
    N0 = Eb./Eb_N0 ;
    k=k+1 ;  
    %##########################
    %AWGN
    noise = (sqrt(N0./2).*randn (t_length,number_of_bits));
    rx_signal = tx_signal+noise ;
    
    %##########################
    %RX
    rx_bit_demodulation = zeros(1,number_of_bits) ;
    for jj = 1:1:number_of_bits
        R1=(sum(rx_signal(:,jj).*cos_11)^2)+(sum(rx_signal(:,jj).*sin_11)^2);
        R2=(sum(rx_signal(:,jj).*cos_22)^2)+(sum(rx_signal(:,jj).*sin_22)^2);
        if (R2>R1)
            rx_bit_demodulation(1,jj) = 1 ;
        end
    end
    %###########################
    %BER
    BER_simulation(1,k) = numel(1,xor(tx_bit_stream,rx_bit_demodulation))./number_of_bits ;
    BER_theorical(1,k) = 0.5*exp(-0.5*Eb_N0);
    BER_theorical_coherent(1,k)=qfunc(sqrt(Eb_N0)) ;
end
%#########################################33
for a = 1:number_of_bits
    for b = 1:t_length
        sent_signal(1,ss)=tx_signal(b,a);
        received_signal(1,ss)=rx_signal(b,a);
        ss=ss+1;
    end
end
%############################
% Ber PLOT
Eb_N0_dB = 0:12 ;
figure(1);
semilogy(Eb_N0_dB,BER_simulation,'b*')
hold on ;
semilogy(Eb_N0_dB,BER_theorical,'r',Eb_N0_dB,BER_theorical_coherent,'k-.')
legend('BER simulation','Noncoherent BER theory','Coherent BER theory') ;
xlabel('Eb/N0 dB') ;ylabel('BER') ;title('Noncoherent FSK') ;grid on ;
%############################
if(t_length*number_of_bits <= 1000)
    figure(2);
    plot(sent_signal);
    xlabel('t*number of bits') ;ylabel('sent signal');
    figure(3);
    plot(received_signal);
    xlabel('t*number of bits') ;ylabel('received signal');
end



