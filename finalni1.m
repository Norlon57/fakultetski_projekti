close all;
clear all;
clc;

s=rng('shuffle');

numFFT = input('Unesite broj FFT odbiraka: '); % Broj FFT odbiraka
m=input('Unesite koliko razlicitih kategorija saobracaja zelite: ');
disp('Unesite koliko blokova resursa zelite da dodelite vasim kategorijama: ');
for i=1:m
    rbSize(1,i)=input('');
end

disp('Unesite broj podnosioca po bloku:');

for i=1:m
    numRBs(1,i)=input('');    % Broj simbola
end

zbir=0;

for i=1:m
    zbir=zbir+rbSize(1,i)*numRBs(1,i);
    if zbir>numFFT
        disp(['Premasili ste opseg od ' num2str(numFFT) ' odbiraka']);
        break;
    end
end

disp('Stavite duzinu ciklicnog prefiksa u %: ');
for i=1:m    
    cpLen(1,i) = input('');         % Duzina ciklicnog prefiksa po simbolu
    cpLen(1,i) = round(numRBs(1,i)/cpLen(1,i)*rbSize(1,i));
end

bitsPerSubCarrier = input('Unesite zeljenu linearnu modulaciju ( 2: QPSK, 4: 16QAM, 6: 64QAM, 8: 256QAM): ');

snrdB = input('Unesite SNR[dB]: ');              % SNR u dB

toneOffset = input('Tonski ofset: ');        % Tonski ofset ili dodatni opseg (u podnosiocima)
L = numFFT/2+1;                % Duzina filtra (=redFiltra+1)

disp(['Iskoristili ste ' num2str(zbir) ' odbiraka od ' num2str(numFFT) ' odbiraka. To je ' num2str(zbir/numFFT*100) ' % ']);

halfFilt = floor(L/2);
n = -halfFilt:halfFilt;

%Sinc truncation prozor
w = (0.5*(1+cos(2*pi.*n/(L-1)))).^0.6;

for i=1:m
    % broj podnosioica u podopsegu
    numDataCarriers=numRBs(1,i).*rbSize(1,i);    
    
    %Sinc funkcija prototip filter
    pb = sinc(((numDataCarriers+2*toneOffset).*n)./numFFT); 
    
    % Koeficijenti normalizovanog NF filtra
    fnum = (pb.*w)/sum(pb).*w;        
    
    % Impulsni odziv filtra
    h = fvtool(fnum, 'Analysis', 'impulse', ...
    'NormalizedFrequency', 'off', 'Fs', 15.36e6);
    h.CurrentAxes.XLabel.String = 'Time (\mus)';
    h.FigureToolbar = 'off';        

    % Upotreba DSP za filtriranje
    filtTx = dsp.FIRFilter('Structure', 'Direct form symmetric', ...
    'Numerator', fnum);
    filtRx = clone(filtTx); % Odgovarajuci filter za prijem (Rx)
    
    % QAM simbolsko mapiranje
    qamMapper = comm.RectangularQAMModulator( ...
    'ModulationOrder', 2^bitsPerSubCarrier, 'BitInput', true, ...
    'NormalizationMethod', 'Average power');

    % Crtanje spektra
    hFig = figure('Position', figposition([46 50 30 30]), 'MenuBar', 'none');
    axis([-0.5 0.5 -200 -20]);
    hold on;
    grid on
    xlabel('Normalizovana frekvencija');
    ylabel('PSD (dBW/Hz)')
    title(['f-OFDM, ' num2str(numRBs(1,i)) ' Podnosioca, '  ...
        num2str(rbSize(1,i)) ' Blokova'])
    legend();
    
    % Generisanje bita
    bitsIn = randi([0 1], bitsPerSubCarrier*numDataCarriers, 1);
    symbolsIn = qamMapper(bitsIn);
    
    % Pakovanje bita u simbole
    offset = (numFFT-numDataCarriers)/2; % za centar opsega
    
    symbolsInOFDM = [zeros(offset,1); symbolsIn; ...
                     zeros(numFFT-offset-numDataCarriers,1)];
    ifftOut = ifft(ifftshift(symbolsInOFDM));
    
    % Dodavanje Ciklicnog prefiksa na pocetak simbola
    txSigOFDM = [ifftOut(end-cpLen(1,i)+1:end); ifftOut];
    
    % Filter, sa dodavanje nula da bi se krajevi odbacili. Uzeti Tx
    txSigFOFDM = filtTx([txSigOFDM; zeros(L-1,1)]);
    
    % Plot spektralna gustina snage (SGS)
    [psd,f] = periodogram(txSigFOFDM, rectwin(length(txSigFOFDM)), ...
                      numFFT*2, 1, 'centered');
    plot(f,10*log10(psd));
        
    % Porednje peak-to-average-power ratio (PAPR)
    PAPR1 = comm.CCDF('PAPROutputPort', true, 'PowerUnits', 'dBW');
    [~,~,paprFOFDM] = PAPR1(txSigFOFDM);
    disp(['PAPR za f-OFDM od ' num2str(numRBs(1,i))  ' = '  num2str(paprFOFDM) ' dB']);

    % Crtanje power spectral density (PSD) za OFDM signal
    [psd,f] = periodogram(txSigOFDM, rectwin(length(txSigOFDM)), numFFT*2, ...
                         1, 'centered');
    hFig = figure('Position', figposition([46 15 30 30]));
    plot(f,10*log10(psd));
    grid on
    axis([-0.5 0.5 -100 -20]);
    xlabel('Normalizovana frekvencija');
    ylabel('PSD (dBW/Hz)')
    title(['OFDM, ' num2str(numRBs(1,i)) ' Podnosioca, '  ...
        num2str(rbSize(1,i)) ' Blokova'])
 
    % Racunaj peak-to-average-power ratio (PAPR)
    PAPR2 = comm.CCDF('PAPROutputPort', true, 'PowerUnits', 'dBW');
    [~,~,paprOFDM] = PAPR2(txSigOFDM);
    disp(['PAPR za OFDM od ' num2str(numRBs(1,i)) ' = ' num2str(paprOFDM) ' dB']);

    % Dodaj aditivni beli gausov sum (AWGN)
    rxSig = awgn(txSigFOFDM, snrdB, 'measured');

    % Prijem odgovarajuceg filtra
    rxSigFilt = filtRx(rxSig);
    
    % Uracunavamo kasnjenje filtra
    rxSigFiltSync = rxSigFilt(L:end);
    
    % Skidanje ciklicnog prefiksa
    rxSymbol = rxSigFiltSync(cpLen(1,i)+1:end);
    
    % FFT
    RxSymbols = fftshift(fft(rxSymbol));
   
    % Odabir podataka podnosioca
    dataRxSymbols = RxSymbols(offset+(1:numDataCarriers));
    
    % Crtanje primljenih simbola u konstelaciju
switch bitsPerSubCarrier
    case 2  % QPSK
        refConst = qammod((0:3).', 4, 'UnitAveragePower', true);
    case 4  % 16QAM
        refConst = qammod((0:15).', 16,'UnitAveragePower', true);
    case 6  % 64QAM
        refConst = qammod((0:63).', 64,'UnitAveragePower', true);
    case 8  % 256QAM
        refConst = qammod((0:255).', 256,'UnitAveragePower', true);
end
    
    % Demapiranje i racunanje BER 
    qamDemod = comm.RectangularQAMDemodulator('ModulationOrder', ...
        2^bitsPerSubCarrier, 'BitOutput', true, ...
       'NormalizationMethod', 'Average power');

    BER = comm.ErrorRate;

    % Izvrsi odairanje i racunaj gresku
    rxBits = qamDemod(dataRxSymbols);
    ber = BER(bitsIn, rxBits);
    
    disp(['f-OFDM Prijem, BER = ' num2str(ber(1)) ' za SNR = ' ...
      num2str(snrdB) ' dB']);
    disp(' ');

end

constDiagRx = comm.ConstellationDiagram( ...
    'ShowReferenceConstellation', true, ...
    'ReferenceConstellation', refConst, ...
    'Position', figposition([20 15 30 40]), ...
    'EnableMeasurements', true, ...
    'MeasurementInterval', length(dataRxSymbols), ...
    'Title', 'F-OFDM Demodulated Symbols', ...
    'Name', 'F-OFDM Reception', ...
    'XLimits', [-1.5 1.5], 'YLimits', [-1.5 1.5]);
constDiagRx(dataRxSymbols) 

rng(s);