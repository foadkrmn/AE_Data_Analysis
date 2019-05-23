% This code is for studying the optimum bin width for histogram based on
% Scott's equation in Scott, D. W. (1979). On optimal and data-based histograms. Biometrika, 66(3), 605-610.
% on valid waveforms used in analysis of test 5RAA07

clear all
clc

% Define file locations
ch1kd = 'G:\My Drive\Research\Fatigue on composites (with Dr Modarres)\TDA Project\Phase 2 Option\Tests\5RAA07_102218\AE\Results\Bin Width Study\Ch1kd.csv';
ch2kd = 'G:\My Drive\Research\Fatigue on composites (with Dr Modarres)\TDA Project\Phase 2 Option\Tests\5RAA07_102218\AE\Results\Bin Width Study\Ch2kd.csv';

skew = 'G:\My Drive\Research\Fatigue on composites (with Dr Modarres)\TDA Project\Phase 2 Option\Tests\5RAA07_102218\AE\Results\Bin Width Study\Skewnes Coefficient.csv';
kurt = 'G:\My Drive\Research\Fatigue on composites (with Dr Modarres)\TDA Project\Phase 2 Option\Tests\5RAA07_102218\AE\Results\Bin Width Study\Kutosis Coefficient.csv';

WaveAddress_root = 'G:\My Drive\Research\Fatigue on composites (with Dr Modarres)\TDA Project\Phase 2 Option\Tests\5RAA07_102218\AE\Waveform\';

% Import data
Ch1_waveforms = csvread(ch1kd);
Ch2_waveforms = csvread(ch2kd);
sk_coeff = csvread(skew);
kur_coeff = csvread(kurt);

% For channel 1 we have:
for k = 1:length(Ch1_waveforms)
    Wave_filename_ch1 = strcat('Test_01_1_',num2str(Ch1_waveforms(k)),'.txt');
    WaveAddress = [WaveAddress_root,Wave_filename_ch1];
    
    
%     waveform file is imported.
    Waveform_DTF_1 = importdata(WaveAddress,' ',12);
    
%    values regarding hit lockout time are all zero. Those values are
%    removed in this step. The criteria for removing is if average value
%    from a specific point to the end is zero, then take the index number
%    and cut the rest out.
       
     for i=1:length(Waveform_DTF_1.data)
         m = mean(Waveform_DTF_1.data(i:length(Waveform_DTF_1.data)));
         if m == 0
             HitLockOutIndex_1(k,:)=i;
             break
         elseif i==10240
             HitLockOutIndex_1(k,:)=i;
             Long_Waveform_ch1(r,:) = Ch1_waveforms(k);
             r=r+1;
         end
     end
     
     if HitLockOutIndex_1(k) ~= length(Waveform_DTF_1.data)
         Waveform_DTF_1.data(HitLockOutIndex_1(k):length(Waveform_DTF_1.data),:)=[];
     end
     
%    Find standard deviation, skewness, kurtosis for the waveform after HLT removal  
     ch1_sd(k,1) = std(Waveform_DTF_1.data);
     ch1_sk(k,1) = skewness(Waveform_DTF_1.data);
     ch1_kur(k,1) = kurtosis(Waveform_DTF_1.data);
     ch1_n(k,1) = length(Waveform_DTF_1.data);
     
%    Calcualte initial bin width before applying correction coefficients  
     ch1_bin(k,1) = 3.49*ch1_sd(k,1)*(ch1_n(k,1).^-0.33);     

end

% Find correction coefficients for each bin width according to its skewness
% and kurtosis
for i=1:length(ch1_bin)
    
%   Find skewness coefficient   
    for j=1:length(sk_coeff)-1
        if abs(ch1_sk(i,1)) < sk_coeff(1,1)
            ch1_sk_coeff(i,1) = sk_coeff(1,2);
        elseif (sk_coeff(j,1) < abs(ch1_sk(i,1)) & abs(ch1_sk(i,1)) < sk_coeff(j+1,1))
            ch1_sk_coeff(i,1) = sk_coeff(j+1,2);
        elseif sk_coeff(length(sk_coeff),1) < abs(ch1_sk(i,1))
            ch1_sk_coeff(i,1) = sk_coeff(length(sk_coeff),2);
        end
    end
    
%   Find kutosis coefficient
    for m=1:length(kur_coeff)-1
        if abs(ch1_kur(i,1)) < kur_coeff(1,1)
            ch1_kur_coeff(i,1) = kur_coeff(1,2);
        elseif (kur_coeff(m,1) < abs(ch1_kur(i,1)) & abs(ch1_kur(i,1)) < kur_coeff(m+1,1))
            ch1_kur_coeff(i,1) = kur_coeff(m+1,2);
        elseif kur_coeff(length(kur_coeff),1) < abs(ch1_kur(i,1))
            ch1_kur_coeff(i,1) = kur_coeff(length(kur_coeff),2);
        end
    end
end

% Calcualte bin width by multiplying coefficients
ch1_bin_mod = ch1_bin.*ch1_sk_coeff.* ch1_kur_coeff;

% Generate a variable containing all usefull informations
ch1_bin_result = [Ch1_waveforms ch1_sd ch1_n ch1_bin ch1_sk ch1_sk_coeff ch1_kur ch1_kur_coeff ch1_bin_mod];


% For channel 2 we have:
for k = 1:length(Ch2_waveforms)
    Wave_filename_ch2 = strcat('Test_01_1_',num2str(Ch2_waveforms(k)),'.txt');
    WaveAddress = [WaveAddress_root,Wave_filename_ch2];
    
    
%     waveform file is imported.
    Waveform_DTF_2 = importdata(WaveAddress,' ',12);
    
%    values regarding hit lockout time are all zero. Those values are
%    removed in this step. The criteria for removing is if average value
%    from a specific point to the end is zero, then take the index number
%    and cut the rest out.
       
     for i=1:length(Waveform_DTF_2.data)
         m = mean(Waveform_DTF_2.data(i:length(Waveform_DTF_2.data)));
         if m == 0
             HitLockOutIndex_2(k,:)=i;
             break
         elseif i==10240
             HitLockOutIndex_2(k,:)=i;
             Long_Waveform_ch2(r,:) = Ch2_waveforms(k);
             r=r+1;
         end
     end
     
     if HitLockOutIndex_2(k) ~= length(Waveform_DTF_2.data)
         Waveform_DTF_2.data(HitLockOutIndex_2(k):length(Waveform_DTF_2.data),:)=[];
     end

%    Find standard deviation, skewness, kurtosis for the waveform after HLT removal  
     ch2_sd(k,1) = std(Waveform_DTF_2.data);
     ch2_sk(k,1) = skewness(Waveform_DTF_2.data);
     ch2_kur(k,1) = kurtosis(Waveform_DTF_2.data);
     ch2_n(k,1) = length(Waveform_DTF_2.data);
     
%    Calcualte initial bin width before applying correction coefficients       
     ch2_bin(k,1) = 3.49*ch2_sd(k,1)*(ch2_n(k,1).^-0.33);     

end

% Find correction coefficients for each bin width according to its skewness
% and kurtosis

for i=1:length(ch2_bin)

    %   Find skewness coefficient   
    for j=1:length(sk_coeff)-1
        if abs(ch2_sk(i,1)) < sk_coeff(1,1)
            ch2_sk_coeff(i,1) = sk_coeff(1,2);
        elseif (sk_coeff(j,1) < abs(ch2_sk(i,1)) & abs(ch2_sk(i,1)) < sk_coeff(j+1,1))
            ch2_sk_coeff(i,1) = sk_coeff(j+1,2);
        elseif sk_coeff(length(sk_coeff),1) < abs(ch2_sk(i,1))
            ch2_sk_coeff(i,1) = sk_coeff(length(sk_coeff),2);
        end
    end
    
    
    %   Find kutosis coefficient
    for m=1:length(kur_coeff)-1
        if abs(ch2_kur(i,1)) < kur_coeff(1,1)
            ch2_kur_coeff(i,1) = kur_coeff(1,2);
        elseif (kur_coeff(m,1) < abs(ch2_kur(i,1)) & abs(ch2_kur(i,1)) < kur_coeff(m+1,1))
            ch2_kur_coeff(i,1) = kur_coeff(m+1,2);
        elseif kur_coeff(length(kur_coeff),1) < abs(ch2_kur(i,1))
            ch2_kur_coeff(i,1) = kur_coeff(length(kur_coeff),2);
        end
    end
end

% Calcualte bin width by multiplying coefficients
ch2_bin_mod = ch2_bin.*ch2_sk_coeff.* ch2_kur_coeff;

% Generate a variable containing all usefull informations
ch2_bin_result = [Ch2_waveforms ch2_sd ch2_n ch2_bin ch2_sk ch2_sk_coeff ch2_kur ch2_kur_coeff ch2_bin_mod];

% Generate excel files with all results
filename1 = 'Ch1 Bin Width Analysis';
filename2 = 'Ch2 Bin Width Analysis';
xlswrite(filename1,ch1_bin_result)
xlswrite(filename2,ch2_bin_result)
