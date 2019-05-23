% % % % % % This part of the code is used to depict signal of same hits in channel 1
% and channel 2 that have the following properties: ent_ch1 > ent_ch2 & is
% the min ent before initiation; ent_ch1 < ent_ch2; ent_ch1 = ent_ch2.
% Waveform and histogram of all three scenarios should be generated.

% scenario 1: ent_ch1 > ent_ch2 & min ent before initiation

S = sort(Ent_dtf_Ch1_2(1:I_crack_init_ch1_AE),'ascend');
for i=1:length(S)
    if S(i)>0
        index1 = find(Ent_dtf_Ch1_2(1:I_crack_init_ch1_AE,1) == S(i));
        break
    end
end

index1_1 = find(Ent_dtf_Ch2_2(:,1) == min(Ent_dtf_Ch2_2(975:1339,1)));

limit_wave_1 = [0 inf -0.1 0.1];
limit_hist_1 = [-0.1 0.1 0 800];

%----------------- second scenario -----------------%
index2(:,1) = find(Ent_dtf_Ch2_2(:,1) > Ent_dtf_Ch1_2(:,1));
index2_selected = index2(811,1);

limit_wave_2 = [0 inf -0.2 0.2];
limit_hist_2 = [-0.2 0.2 0 800];

%----------------- third scenario -----------------%
index3(:,1) = find(abs(Ent_dtf_Ch2_2(:,1) - Ent_dtf_Ch1_2(:,1)) <= 0.001);
index3_selected = index3(4,1);

limit_wave_3 = [0 inf -0.15 0.15];
limit_hist_3 = [-0.15 0.15 0 600];



r=1;
% For channel 1 we have:
for k = index1_1;
    
    Wave_filename_ch1 = strcat('Test_01_1_',num2str(ch1kd(k)),'.txt');
    WaveAddress = [WaveAddress_root,Wave_filename_ch1];

    
% % % % % % % % % % This part of the code is the one reffered to as updated code % % % % % % % % % %       
    % First check whether the hit waveform is selected right, by comparing
    % arrival time on waveform and in ch1dtf(:,2)
    Waveform_DTF_1_t = importdata(WaveAddress,' ',10);
    
    y=1;
    if Waveform_DTF_1_t.data < ch1dtf(k,2)
        while Waveform_DTF_1_t.data < ch1dtf(k,2)
            Wave_filename_ch1 = strcat('Test_01_1_',num2str(ch1kd(k)+y),'.txt');
            WaveAddress = [WaveAddress_root,Wave_filename_ch1];
            Waveform_DTF_1_t = importdata(WaveAddress,' ',10);
            y=y+1;
        end
        disp(['Waveform Test_01_1_',num2str(ch1kd(k)+y),' is used instead of',' Test_01_1_',num2str(ch1kd(k)),'.'])
    end 
    
    y=1;
    if Waveform_DTF_1_t.data > ch1dtf(k,2)
        while Waveform_DTF_1_t.data > ch1dtf(k,2)
            Wave_filename_ch1 = strcat('Test_01_1_',num2str(ch1kd(k)-y),'.txt');
            WaveAddress = [WaveAddress_root,Wave_filename_ch1];
            Waveform_DTF_1_t = importdata(WaveAddress,' ',10);
            y=y+1;
        end
        disp(['Waveform Test_01_1_',num2str(ch1kd(k)-y),' is used instead of',' Test_01_1_',num2str(ch1kd(k)),'.'])
    end 
    
    if Waveform_DTF_1_t.data ~= ch1dtf(k,2)
        disp(['No waveform found for hit index number ',num2str(k),' in ch2dtf data.'])
        break
    end
% % % % % % % % % % This part of the code is the one reffered to as updated code % % % % % % % % % %       
  

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
             Long_Waveform_ch1(r,:) = ch1kd(k);
             r=r+1;
         end
     end
     
     if HitLockOutIndex_1(k) ~= length(Waveform_DTF_1.data)
         Waveform_DTF_1.data(HitLockOutIndex_1(k):length(Waveform_DTF_1.data),:)=[];
     end
      
%     Plot waveform
      f25 = figure('units','normalized','outerposition',[0 0 0.6 0.7]);
      plot(Waveform_DTF_1.data)
      title (['Ch1 HLT removed Signal at cycle number: ',num2str(Ent_dtf_Ch1_2(k,2))])
      ylabel('Voltage(volts)')
      axis(limit_wave_1)
      legend(['Signal waveform after HTL removal, ','Signal Amplitude is: ',num2str(ch1dtf(k,10))])
      print(f25,'01_sig_01_ch1_waveform.tiff','-dtiffn','-r300')         
     

      f22 = figure('units','normalized','outerposition',[0 0 0.6 0.7]);
      histogram(Waveform_DTF_1.data,edge_2)
      xlabel('Voltage(volts)')
      axis(limit_hist_1)
      title(['Waveform histogram - Ch1 - ','cycle number: ',num2str(Ent_dtf_Ch1_2(k,2))])
      legend(['bin width = ',num2str(BinWidth2),'  ','Entropy is: ',num2str(Ent_dtf_Ch1_2(k,1))])
      print(f22,'03_sig_01_ch1_hist.tiff','-dtiffn','-r300')
      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Channel 2 starts here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    Wave_filename_ch2 = strcat('Test_01_2_',num2str(ch2kd(k)),'.txt');
    WaveAddress = [WaveAddress_root,Wave_filename_ch2];
    
% % % % % % % % % % This part of the code is the one reffered to as updated code % % % % % % % % % %    
    % First check whether the hit waveform is selected right, by comparing
    % arrival time on waveform and in ch1dtf(:,2)
    Waveform_DTF_2_t = importdata(WaveAddress,' ',10);
    
    y=1;
    if Waveform_DTF_2_t.data < ch2dtf(k,2)
        while Waveform_DTF_2_t.data < ch2dtf(k,2)
            Wave_filename_ch2 = strcat('Test_01_2_',num2str(ch2kd(k)+y),'.txt');
            WaveAddress = [WaveAddress_root,Wave_filename_ch2];
            Waveform_DTF_2_t = importdata(WaveAddress,' ',10);
            y=y+1;
        end
        disp(['Waveform Test_01_2_',num2str(ch2kd(k)+y),' is used instead of',' Test_01_2_',num2str(ch2kd(k)),'.'])
    end 
    
    y=1;
    if Waveform_DTF_2_t.data > ch2dtf(k,2)
        while Waveform_DTF_2_t.data > ch2dtf(k,2)
            Wave_filename_ch2 = strcat('Test_01_2_',num2str(ch2kd(k)-y),'.txt');
            WaveAddress = [WaveAddress_root,Wave_filename_ch2];
            Waveform_DTF_2_t = importdata(WaveAddress,' ',10);
            y=y+1;
        end
        disp(['Waveform Test_01_2_',num2str(ch2kd(k)-y),' is used instead of',' Test_01_2_',num2str(ch2kd(k)),'.'])
    end 
    
    if Waveform_DTF_2_t.data ~= ch2dtf(k,2)
        disp(['No waveform found for hit index number ',num2str(k),' in ch2dtf data.'])
        break
    end
% % % % % % % % % % This part of the code is the one reffered to as updated code % % % % % % % % % %       
  

%     waveform file is imported.
    Waveform_DTF_2 = importdata(WaveAddress,' ',12);
    
%    values regarding hit lockout time are all zero. Those values are
%    removed in this step. The criteria for removing is if average value
%    from a specific point to the end is zero, then take the index number
%    and cut the rest out.
     for i=1:length(Waveform_DTF_2.data)
         m = mean(Waveform_DTF_2.data(i:length(Waveform_DTF_2.data)));
         if m == 0
             HitLockOutIndex_1(k,:)=i;
             break
         elseif i==10240
             HitLockOutIndex_1(k,:)=i;
             Long_Waveform_ch2(r,:) = ch2kd(k);
             r=r+1;
         end
     end
     
     if HitLockOutIndex_1(k) ~= length(Waveform_DTF_2.data)
         Waveform_DTF_2.data(HitLockOutIndex_1(k):length(Waveform_DTF_2.data),:)=[];
     end
      
%     Plot waveform
      f25 = figure('units','normalized','outerposition',[0 0 0.6 0.7]);
      plot(Waveform_DTF_2.data)
      title (['Ch2 HLT removed Signal at cycle number: ',num2str(Ent_dtf_Ch2_2(k,2))])
      ylabel('Voltage(volts)')
      axis(limit_wave_1)
      legend(['Signal waveform after HTL removal, ','Signal Amplitude is: ',num2str(ch2dtf(k,10))])
      print(f25,'02_sig_01_ch2_waveform.tiff','-dtiffn','-r300')         
     
      f22 = figure('units','normalized','outerposition',[0 0 0.6 0.7]);
      histogram(Waveform_DTF_2.data,edge_2)
      xlabel('Voltage(volts)')
      axis(limit_hist_1)
      title(['Waveform histogram - Ch2 - ','cycle number: ',num2str(Ent_dtf_Ch2_2(k,2))])
      legend(['bin width = ',num2str(BinWidth2),'  ','Entropy is: ',num2str(Ent_dtf_Ch2_2(k,1))])
      print(f22,'04_sig_01_ch2_hist.tiff','-dtiffn','-r300')
    
end









%----------------------% scenario 2: ent_ch1 < ent_ch2 %-------------------------------



r=1;
% For channel 1 we have:
for k = index2_selected;
    
    Wave_filename_ch1 = strcat('Test_01_1_',num2str(ch1kd(k)),'.txt');
    WaveAddress = [WaveAddress_root,Wave_filename_ch1];

    
% % % % % % % % % % This part of the code is the one reffered to as updated code % % % % % % % % % %       
    % First check whether the hit waveform is selected right, by comparing
    % arrival time on waveform and in ch1dtf(:,2)
    Waveform_DTF_1_t = importdata(WaveAddress,' ',10);
    
    y=1;
    if Waveform_DTF_1_t.data < ch1dtf(k,2)
        while Waveform_DTF_1_t.data < ch1dtf(k,2)
            Wave_filename_ch1 = strcat('Test_01_1_',num2str(ch1kd(k)+y),'.txt');
            WaveAddress = [WaveAddress_root,Wave_filename_ch1];
            Waveform_DTF_1_t = importdata(WaveAddress,' ',10);
            y=y+1;
        end
        disp(['Waveform Test_01_1_',num2str(ch1kd(k)+y),' is used instead of',' Test_01_1_',num2str(ch1kd(k)),'.'])
    end 
    
    y=1;
    if Waveform_DTF_1_t.data > ch1dtf(k,2)
        while Waveform_DTF_1_t.data > ch1dtf(k,2)
            Wave_filename_ch1 = strcat('Test_01_1_',num2str(ch1kd(k)-y),'.txt');
            WaveAddress = [WaveAddress_root,Wave_filename_ch1];
            Waveform_DTF_1_t = importdata(WaveAddress,' ',10);
            y=y+1;
        end
        disp(['Waveform Test_01_1_',num2str(ch1kd(k)-y),' is used instead of',' Test_01_1_',num2str(ch1kd(k)),'.'])
    end 
    
    if Waveform_DTF_1_t.data ~= ch1dtf(k,2)
        disp(['No waveform found for hit index number ',num2str(k),' in ch2dtf data.'])
        break
    end
% % % % % % % % % % This part of the code is the one reffered to as updated code % % % % % % % % % %       
  

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
             Long_Waveform_ch1(r,:) = ch1kd(k);
             r=r+1;
         end
     end
     
     if HitLockOutIndex_1(k) ~= length(Waveform_DTF_1.data)
         Waveform_DTF_1.data(HitLockOutIndex_1(k):length(Waveform_DTF_1.data),:)=[];
     end
      
%     Plot waveform
      f25 = figure('units','normalized','outerposition',[0 0 0.6 0.7]);
      plot(Waveform_DTF_1.data)
      title (['Ch1 HLT removed Signal at cycle number: ',num2str(ch1dtf(find(ch1dtf(:,18) == ch1kd(k)),19))])
      ylabel('Voltage(volts)')
      axis(limit_wave_2)
      legend(['Signal waveform after HTL removal, ','Signal Amplitude is: ',num2str(ch1dtf(k,10))])
      print(f25,'05_sig_02_ch1_waveform.tiff','-dtiffn','-r300')         
     
      f22 = figure('units','normalized','outerposition',[0 0 0.6 0.7]);
      histogram(Waveform_DTF_1.data,edge_2)
      xlabel('Voltage(volts)')
      axis(limit_hist_2)
      title(['Waveform histogram - Ch1 - ','cycle number: ',num2str(ch1dtf(find(ch1dtf(:,18) == ch1kd(k)),19))])
      legend(['bin width = ',num2str(BinWidth2),'  ','Entropy is: ',num2str(Ent_dtf_Ch1_2(k,1))])
      print(f22,'07_sig_02_ch1_hist.tiff','-dtiffn','-r300')

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Channel 2 starts here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
Wave_filename_ch2 = strcat('Test_01_2_',num2str(ch2kd(k)),'.txt');
    WaveAddress = [WaveAddress_root,Wave_filename_ch2];
    
% % % % % % % % % % This part of the code is the one reffered to as updated code % % % % % % % % % %    
    % First check whether the hit waveform is selected right, by comparing
    % arrival time on waveform and in ch1dtf(:,2)
    Waveform_DTF_2_t = importdata(WaveAddress,' ',10);
    
    y=1;
    if Waveform_DTF_2_t.data < ch2dtf(k,2)
        while Waveform_DTF_2_t.data < ch2dtf(k,2)
            Wave_filename_ch2 = strcat('Test_01_2_',num2str(ch2kd(k)+y),'.txt');
            WaveAddress = [WaveAddress_root,Wave_filename_ch2];
            Waveform_DTF_2_t = importdata(WaveAddress,' ',10);
            y=y+1;
        end
        disp(['Waveform Test_01_2_',num2str(ch2kd(k)+y),' is used instead of',' Test_01_2_',num2str(ch2kd(k)),'.'])
    end 
    
    y=1;
    if Waveform_DTF_2_t.data > ch2dtf(k,2)
        while Waveform_DTF_2_t.data > ch2dtf(k,2)
            Wave_filename_ch2 = strcat('Test_01_2_',num2str(ch2kd(k)-y),'.txt');
            WaveAddress = [WaveAddress_root,Wave_filename_ch2];
            Waveform_DTF_2_t = importdata(WaveAddress,' ',10);
            y=y+1;
        end
        disp(['Waveform Test_01_2_',num2str(ch2kd(k)-y),' is used instead of',' Test_01_2_',num2str(ch2kd(k)),'.'])
    end 
    
    if Waveform_DTF_2_t.data ~= ch2dtf(k,2)
        disp(['No waveform found for hit index number ',num2str(k),' in ch2dtf data.'])
        break
    end
% % % % % % % % % % This part of the code is the one reffered to as updated code % % % % % % % % % %       
  

%     waveform file is imported.
    Waveform_DTF_2 = importdata(WaveAddress,' ',12);
    
%    values regarding hit lockout time are all zero. Those values are
%    removed in this step. The criteria for removing is if average value
%    from a specific point to the end is zero, then take the index number
%    and cut the rest out.
     for i=1:length(Waveform_DTF_2.data)
         m = mean(Waveform_DTF_2.data(i:length(Waveform_DTF_2.data)));
         if m == 0
             HitLockOutIndex_1(k,:)=i;
             break
         elseif i==10240
             HitLockOutIndex_1(k,:)=i;
             Long_Waveform_ch2(r,:) = ch2kd(k);
             r=r+1;
         end
     end
     
     if HitLockOutIndex_1(k) ~= length(Waveform_DTF_2.data)
         Waveform_DTF_2.data(HitLockOutIndex_1(k):length(Waveform_DTF_2.data),:)=[];
     end
      
%     Plot waveform
      f25 = figure('units','normalized','outerposition',[0 0 0.6 0.7]);
      plot(Waveform_DTF_2.data)
      title (['Ch2 HLT removed Signal at cycle number: ',num2str(ch2dtf(find(ch2dtf(:,18) == ch2kd(k)),19))])
      ylabel('Voltage(volts)')
      axis(limit_wave_2)
      legend(['Signal waveform after HTL removal, ','Signal Amplitude is: ',num2str(ch2dtf(k,10))])
      print(f25,'06_sig_02_ch2_waveform.tiff','-dtiffn','-r300')         
     
      f22 = figure('units','normalized','outerposition',[0 0 0.6 0.7]);
      histogram(Waveform_DTF_2.data,edge_2)
      xlabel('Voltage(volts)')
      axis(limit_hist_2)
      title(['Waveform histogram - Ch2 - ','cycle number: ',num2str(ch2dtf(find(ch2dtf(:,18) == ch2kd(k)),19))])
      legend(['bin width = ',num2str(BinWidth2),'  ','Entropy is: ',num2str(Ent_dtf_Ch2_2(k,1))])
      print(f22,'08_sig_02_ch2_hist.tiff','-dtiffn','-r300')
    
end










% --------------------------% scenario 3: ent_ch1 = ent_ch2 %-----------------------------


r=1;
% For channel 1 we have:
for k = index3_selected;

    Wave_filename_ch1 = strcat('Test_01_1_',num2str(ch1kd(k)),'.txt');
    WaveAddress = [WaveAddress_root,Wave_filename_ch1];

    
% % % % % % % % % % This part of the code is the one reffered to as updated code % % % % % % % % % %       
    % First check whether the hit waveform is selected right, by comparing
    % arrival time on waveform and in ch1dtf(:,2)
    Waveform_DTF_1_t = importdata(WaveAddress,' ',10);
    
    y=1;
    if Waveform_DTF_1_t.data < ch1dtf(k,2)
        while Waveform_DTF_1_t.data < ch1dtf(k,2)
            Wave_filename_ch1 = strcat('Test_01_1_',num2str(ch1kd(k)+y),'.txt');
            WaveAddress = [WaveAddress_root,Wave_filename_ch1];
            Waveform_DTF_1_t = importdata(WaveAddress,' ',10);
            y=y+1;
        end
        disp(['Waveform Test_01_1_',num2str(ch1kd(k)+y),' is used instead of',' Test_01_1_',num2str(ch1kd(k)),'.'])
    end 
    
    y=1;
    if Waveform_DTF_1_t.data > ch1dtf(k,2)
        while Waveform_DTF_1_t.data > ch1dtf(k,2)
            Wave_filename_ch1 = strcat('Test_01_1_',num2str(ch1kd(k)-y),'.txt');
            WaveAddress = [WaveAddress_root,Wave_filename_ch1];
            Waveform_DTF_1_t = importdata(WaveAddress,' ',10);
            y=y+1;
        end
        disp(['Waveform Test_01_1_',num2str(ch1kd(k)-y),' is used instead of',' Test_01_1_',num2str(ch1kd(k)),'.'])
    end 
    
    if Waveform_DTF_1_t.data ~= ch1dtf(k,2)
        disp(['No waveform found for hit index number ',num2str(k),' in ch2dtf data.'])
        break
    end
% % % % % % % % % % This part of the code is the one reffered to as updated code % % % % % % % % % %      
  

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
             Long_Waveform_ch1(r,:) = ch1kd(k);
             r=r+1;
         end
     end
     
     if HitLockOutIndex_1(k) ~= length(Waveform_DTF_1.data)
         Waveform_DTF_1.data(HitLockOutIndex_1(k):length(Waveform_DTF_1.data),:)=[];
     end
      
%     Plot waveform
      f25 = figure('units','normalized','outerposition',[0 0 0.6 0.7]);
      plot(Waveform_DTF_1.data)
      title (['Ch1 HLT removed Signal at cycle number: ',num2str(ch1dtf(find(ch1dtf(:,18) == ch1kd(k)),19))])
      ylabel('Voltage(volts)')
      axis(limit_wave_3)
      legend(['Signal waveform after HTL removal, ','Signal Amplitude is: ',num2str(ch1dtf(k,10))])
      print(f25,'09_sig_03_ch1_waveform.tiff','-dtiffn','-r300')         
     
    f22 = figure('units','normalized','outerposition',[0 0 0.6 0.7]);
    histogram(Waveform_DTF_1.data,edge_2)
    xlabel('Voltage(volts)')
    axis(limit_hist_3)
    title(['Waveform histogram - Ch1 - ','cycle number: ',num2str(ch1dtf(find(ch1dtf(:,18) == ch1kd(k)),19))])
    legend(['bin width = ',num2str(BinWidth2),'  ','Entropy is: ',num2str(Ent_dtf_Ch1_2(k,1))])
    print(f22,'11_sig_03_ch1_hist.tiff','-dtiffn','-r300')

    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Channel 2 starts here %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    
    Wave_filename_ch2 = strcat('Test_01_2_',num2str(ch2kd(k)),'.txt');
    WaveAddress = [WaveAddress_root,Wave_filename_ch2];
    
% % % % % % % % % % This part of the code is the one reffered to as updated code % % % % % % % % % %    
    % First check whether the hit waveform is selected right, by comparing
    % arrival time on waveform and in ch1dtf(:,2)
    Waveform_DTF_2_t = importdata(WaveAddress,' ',10);
    
    y=1;
    if Waveform_DTF_2_t.data < ch2dtf(k,2)
        while Waveform_DTF_2_t.data < ch2dtf(k,2)
            Wave_filename_ch2 = strcat('Test_01_2_',num2str(ch2kd(k)+y),'.txt');
            WaveAddress = [WaveAddress_root,Wave_filename_ch2];
            Waveform_DTF_2_t = importdata(WaveAddress,' ',10);
            y=y+1;
        end
        disp(['Waveform Test_01_2_',num2str(ch2kd(k)+y),' is used instead of',' Test_01_2_',num2str(ch2kd(k)),'.'])
    end 
    
    y=1;
    if Waveform_DTF_2_t.data > ch2dtf(k,2)
        while Waveform_DTF_2_t.data > ch2dtf(k,2)
            Wave_filename_ch2 = strcat('Test_01_2_',num2str(ch2kd(k)-y),'.txt');
            WaveAddress = [WaveAddress_root,Wave_filename_ch2];
            Waveform_DTF_2_t = importdata(WaveAddress,' ',10);
            y=y+1;
        end
        disp(['Waveform Test_01_2_',num2str(ch2kd(k)-y),' is used instead of',' Test_01_2_',num2str(ch2kd(k)),'.'])
    end 
    
    if Waveform_DTF_2_t.data ~= ch2dtf(k,2)
        disp(['No waveform found for hit index number ',num2str(k),' in ch2dtf data.'])
        break
    end
% % % % % % % % % % This part of the code is the one reffered to as updated code % % % % % % % % % %       
  

%     waveform file is imported.
    Waveform_DTF_2 = importdata(WaveAddress,' ',12);
    
%    values regarding hit lockout time are all zero. Those values are
%    removed in this step. The criteria for removing is if average value
%    from a specific point to the end is zero, then take the index number
%    and cut the rest out.
     for i=1:length(Waveform_DTF_2.data)
         m = mean(Waveform_DTF_2.data(i:length(Waveform_DTF_2.data)));
         if m == 0
             HitLockOutIndex_1(k,:)=i;
             break
         elseif i==10240
             HitLockOutIndex_1(k,:)=i;
             Long_Waveform_ch2(r,:) = ch2kd(k);
             r=r+1;
         end
     end
     
     if HitLockOutIndex_1(k) ~= length(Waveform_DTF_2.data)
         Waveform_DTF_2.data(HitLockOutIndex_1(k):length(Waveform_DTF_2.data),:)=[];
     end
      
%     Plot waveform
      f25 = figure('units','normalized','outerposition',[0 0 0.6 0.7]);
      plot(Waveform_DTF_2.data)
      title (['Ch2 HLT removed Signal at cycle number: ',num2str(ch2dtf(find(ch2dtf(:,18) == ch2kd(k)),19))])
      ylabel('Voltage(volts)')
      axis(limit_wave_3)
      legend(['Signal waveform after HTL removal, ','Signal Amplitude is: ',num2str(ch2dtf(k,10))])
      print(f25,'10_sig_03_ch2_waveform.tiff','-dtiffn','-r300')         
     

      f22 = figure('units','normalized','outerposition',[0 0 0.6 0.7]);
      histogram(Waveform_DTF_2.data,edge_2)
      xlabel('Voltage(volts)')
      axis(limit_hist_3)
      title(['Waveform histogram - Ch2 - ','cycle number: ',num2str(ch2dtf(find(ch2dtf(:,18) == ch2kd(k)),19))])
      legend(['bin width = ',num2str(BinWidth2),'  ','Entropy is: ',num2str(Ent_dtf_Ch2_2(k,1))])
      print(f22,'12_sig_03_ch2_hist.tiff','-dtiffn','-r300')
    
end
