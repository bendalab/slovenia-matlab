%Autor Matthias Ph. Baumann
%Dieses Skript liest die onset daten die in einem cell angeordnet sein
%müssen, wobei jede 'Spalte' einem Trace entspricht,ein. führt eine
%kreuzkorrelation durch und plottet das Histogramm der zeitdifferenzen.
%Eine weitere Kreuzcorrelatin (mit xcorr) wird mit einer aus den Daten
%generierten time series die mit einem gauss gefaltet ist durchgeführt. 
%Mittels xcov wir der correlationscoeffizient ermittelt.

clear
close all
% Read in Onset-Files (hopefully all files in current folder)
mat = dir('*.mat');
for q = 1:length(mat)
    data{q}=load(mat(q).name);
end

% Prepare space for crosscorrelation coeficients
rxcoef12all=nan(length(mat),1);
rxcoef23all=nan(length(mat),1);
rxcoef13all=nan(length(mat),1);
% main loop over all recordings
for ii=1:length(mat)
    %VERY IMPORTANT clear old traces cause they can differ in size for each
    %loop
    clear datatrace1 datatrace2 datatrace3 corrdata
    % load the 3 traces (spcified for how i saved my onset data)
    load(mat(ii,1).name)
    datatrace1=onset{1,1};
    datatrace2=onset{1,2};
    datatrace3=onset{1,3};
    datatraces={datatrace1,datatrace2,datatrace3};
    % Manual cross correlation
    for k=1:3
        if k==3
            for kk=1:length(datatraces{k-2})
                for kkk=1:length(datatraces{k})
                    corrdata{k}(kk,kkk)=datatraces{k-2}(kk)-datatraces{k}(kkk);
                end
            end
        else
            for kk=1:length(datatraces{k})
                for kkk=1:length(datatraces{k+1})
                    corrdata{k}(kk,kkk)=datatraces{k}(kk)-datatraces{k+1}(kkk);
                end
            end
        end
    end
    % Plot
    figure
    binwidthh=1;%seconds
    for iii=1:length(corrdata)
        subplot(3,1,iii)
        hold on
        histogram(corrdata{iii},'binwidth',binwidthh)
        xlim([-200 200])
        ylim([0 30])
        hold off
    end
    title(mat(ii,1).name)%name of the file 
    
    % I reduced the sample size from 50000 Hz to 500Hz to make it
    % calculateable for matlab in a reasonable time
    samplerate=1/500;
    hertz=500;
    % clear variables from previous Loop.Very important again
    clear maxall zeitraum trace1logic trace2logic trace3logic
    maxall=max([max(datatrace1),max(datatrace2),max(datatrace3)]);
    zeitraum=0:samplerate:round(maxall);
    trace1logic=zeros(length(zeitraum),1);
    trace2logic=zeros(length(zeitraum),1);
    trace3logic=zeros(length(zeitraum),1);
    trace1logic(round(datatrace1.*hertz))=1;
    trace2logic(round(datatrace2.*hertz))=1;
    trace3logic(round(datatrace3.*hertz))=1;
    % Gaussian for Convolution
    sigma=0.5;
    width=4;
    tgauss=-width*sigma:samplerate:width*sigma;
    gaussfilter=exp(-0.5*(tgauss./sigma).^2) *1/(sigma*sqrt(2*pi));
    %Again very very Important
    clear trace1gefaltet trace3gefaltet trace2gefaltet rcc12 rcc23 rcc13 lagscc12 lagscc23 lagscc13 ...
        rxcoef12mat rxcoef12 rxcoef23mat rxcoef23 rxcoef13mat rxcoef13 rxcoef12idx rxcoef23idx rxcoef13idx...
        lengthdiff12 lengthdiff23 lengthdiff13
    %Convolution
    trace1gefaltet=conv(trace1logic,gaussfilter);
    trace2gefaltet=conv(trace2logic,gaussfilter);
    trace3gefaltet=conv(trace3logic,gaussfilter);
    %Calculate the Correlationcoefficent for trace 1 and 2, 2 and 3, 1 and
    %3 I used xcov(Calculates crosscovarianz and substracts the means) to
    %calculate it. The results look reasonable but it is strange how the
    %normalisation with the standartdevariation looks so i did not do this
    if length(trace1gefaltet)==length(trace2gefaltet)
        [rcc12,lagscc12] = xcorr(trace1gefaltet,trace2gefaltet,'coeff');% coeff means that the hight is normalized--> the value at zero from the autocorrelation is set to 1
        rxcoef12mat=xcov(trace1gefaltet,trace2gefaltet,'coeff');
        rxcoef12idx=find(rxcoef12mat==(max(abs(rxcoef12mat))));
        rxcoef12=rxcoef12mat(rxcoef12idx);%/ (std(trace1gefaltet)*std(trace2gefaltet));
    elseif length(trace1gefaltet)<length(trace2gefaltet)
        lengthdiff12=length(trace2gefaltet)-length(trace1gefaltet);
        trace1gefaltet=[trace1gefaltet;zeros(lengthdiff12,1)];
        [rcc12,lagscc12] = xcorr(trace1gefaltet,trace2gefaltet,'coeff');
        rxcoef12mat=xcov(trace1gefaltet,trace2gefaltet,'coeff');
        rxcoef12idx=find(rxcoef12mat==(max(abs(rxcoef12mat))));
        rxcoef12=rxcoef12mat(rxcoef12idx);%/ (std(trace1gefaltet)*std(trace2gefaltet(1:length(trace1gefaltet))));
    else
        lengthdiff12=length(trace1gefaltet)-length(trace2gefaltet);
        trace2gefaltet=[trace2gefaltet;zeros(lengthdiff12,1)];
        [rcc12,lagscc12] = xcorr(trace1gefaltet,trace2gefaltet,'coeff');
        rxcoef12mat=xcov(trace1gefaltet,trace2gefaltet,'coeff');
        rxcoef12idx=find(rxcoef12mat==(max(abs(rxcoef12mat))));
        rxcoef12=rxcoef12mat(rxcoef12idx);%/ (std(trace1gefaltet(1:length(trace2gefaltet)))*std(trace2gefaltet));
    end
    
    if length(trace2gefaltet)==length(trace3gefaltet)
        [rcc23,lagscc23] = xcorr(trace2gefaltet,trace3gefaltet,'coeff');
        rxcoef23mat=xcov(trace2gefaltet,trace3gefaltet,'coeff');
        rxcoef23idx=find(rxcoef23mat==(max(abs(rxcoef23mat))));
        rxcoef23=rxcoef23mat(rxcoef23idx);%/ (std(trace2gefaltet)*std(trace3gefaltet));
        
    elseif length(trace2gefaltet)<length(trace3gefaltet)
        lengthdiff23=length(trace3gefaltet)-length(trace2gefaltet);
        trace2gefaltet=[trace2gefaltet;zeros(lengthdiff23,1)];
        [rcc23,lagscc23] = xcorr(trace2gefaltet,trace3gefaltet,'coeff');
        rxcoef23mat=xcov(trace2gefaltet,trace3gefaltet,'coeff');
        rxcoef23idx=find(rxcoef23mat==(max(abs(rxcoef23mat))));
        rxcoef23=rxcoef23mat(rxcoef23idx);%/ (std(trace2gefaltet)*std(trace3gefaltet(1:length(trace2gefaltet))));
    else
        lengthdiff23=length(trace2gefaltet)-length(trace3gefaltet);
        trace3gefaltet=[trace3gefaltet;zeros(lengthdiff23,1)];
        [rcc23,lagscc23] = xcorr(trace2gefaltet,trace3gefaltet,'coeff');
        rxcoef23mat=xcov(trace2gefaltet,trace3gefaltet,'coeff');
        rxcoef23idx=find(rxcoef23mat==(max(abs(rxcoef23mat))));
        rxcoef23=rxcoef23mat(rxcoef23idx);%/ (std(trace2gefaltet(1:length(trace3gefaltet)))*std(trace3gefaltet));
    end
    
    if length(trace1gefaltet)==length(trace3gefaltet)
        [rcc13,lagscc13] = xcorr(trace1gefaltet,trace3gefaltet,'coeff');
        rxcoef13mat=xcov(trace1gefaltet,trace3gefaltet,'coeff');
        rxcoef13idx=find(rxcoef13mat==(max(abs(rxcoef13mat))));
        rxcoef13=rxcoef13mat(rxcoef13idx);%/ (std(trace1gefaltet)*std(trace3gefaltet));
    elseif length(trace1gefaltet)<length(trace3gefaltet)
        lengthdiff13=length(trace3gefaltet)-length(trace1gefaltet);
        trace1gefaltet=[trace1gefaltet;zeros(lengthdiff13,1)];
        [rcc13,lagscc13]= xcorr(trace1gefaltet,trace3gefaltet,'coeff');
        rxcoef13mat=xcov(trace1gefaltet,trace3gefaltet,'coeff');
        rxcoef13idx=find(rxcoef13mat==(max(abs(rxcoef13mat))));
        rxcoef13=rxcoef13mat(rxcoef13idx);%/ (std(trace1gefaltet)*std(trace3gefaltet(1:length(trace1gefaltet))));
    else
        lengthdiff13=length(trace1gefaltet)-length(trace3gefaltet);
        trace3gefaltet=[trace3gefaltet;zeros(lengthdiff13,1)];
        [rcc13,lagscc13] = xcorr(trace1gefaltet,trace3gefaltet,'coeff');
        rxcoef13mat=xcov(trace1gefaltet,trace3gefaltet,'coeff');
        rxcoef13idx=find(rxcoef13mat==(max(abs(rxcoef13mat))));
        rxcoef13=rxcoef13mat(rxcoef13idx);%/ (std(trace1gefaltet(1:length(trace3gefaltet)))*std(trace3gefaltet));
    end
    
    % plot again the results of the crosscorrelation from the convoluted
    % data with the coefficents in the title the upper plot is always trace
    % 1and2 than 2and3 than 1and3
    
    figure
    subplot(3,1,1)
    plot(lagscc12/hertz,rcc12,'r')
    titl12=rxcoef12;
    title(titl12)
    subplot(3,1,2)
    plot(lagscc23/hertz,rcc23,'g')
    titl23=rxcoef23;
    
    title(titl23)
    subplot(3,1,3)
    plot(lagscc13/hertz,rcc13,'b')
    titl13=rxcoef13;
    
    title(titl13)
    
    % collect the coeficents
    if isempty(rxcoef12)==1
        rxcoef12all(ii)=nan;
    else
        rxcoef12all(ii)=rxcoef12;
    end
    if isempty(rxcoef23)==1
        rxcoef23all(ii)=nan;
    else
        rxcoef23all(ii)=rxcoef23;
    end
    if isempty(rxcoef13)==1
        rxcoef13all(ii)=nan;
    else
        rxcoef13all(ii)=rxcoef13;
    end
    
end