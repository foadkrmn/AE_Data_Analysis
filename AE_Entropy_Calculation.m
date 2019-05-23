% This code calculates waveform entropy for the hits selected using
% AE_Entropy_hit_filter code
clear all; clc;

disp(datetime('now'));
tic;

% Specify all file locations here for convenience. 
% Address_root = 'D:\TDA Project\Phase 2 Option\Tests\CTAA04_020419\AE\';
% Address_root = 'D:\TDA Project\Phase 2 Option\Tests\CTAA05_020519\AE\';
% Address_root = 'D:\TDA Project\Phase 2 Option\Tests\CTAA06_020619\AE\';
% Address_root = 'D:\TDA Project\Phase 2 Option\Tests\CTAA07_020619\AE\';
% Address_root = 'C:\Users\foadkrmn\Documents\Foad\TDA\CTAA08_020719\AE\';
% Address_root = 'C:\Users\foadkrmn\Documents\Foad\TDA\CTAA09_021919\AE\';
% Address_root = 'D:\Foad\TDA\CTAA10_022019\AE\';
% Address_root = 'D:\Foad\TDA\CTAA12_022219\AE\';
% Address_root = 'D:\Foad\TDA\CTAA13_022519\AE\';
% Address_root = 'D:\Foad\TDA\CTAA14_022619\AE\';
Address_root = 'C:\Users\foadkrmn\Documents\Foad\TDA\CTAA15_022719\AE\';
% Address_root = 'C:\Users\foadkrmn\Documents\Foad\TDA\CTAA16_022819\AE\';

% Specify number of Para data files for each test.
para_num = 12; % CTAA04: 10 -- CTAA05: 1  -- CTAA06: 2  -- CTAA07: 2  -- CTAA08: 8
               % CTAA09: 11 -- CTAA10: 13 -- CTAA11: 15 -- CTAA12: 11
               % CTAA13: 11 -- CTAA14: 13 -- CTAA15: 12 -- CTAA16: 11

% Import hit filtered data
HitData_Ch1 = importdata([Address_root,'Results\Ch1_load_deltaT_filtered.txt'],' ',1);
HitData_Ch2 = importdata([Address_root,'Results\Ch2_load_deltaT_filtered.txt'],' ',1);

% Import Sum_growth_## files
for i=1:para_num
    
    if i < 10
        sum_filename = ['Sum_growth_0',num2str(i),'.TXT'];
    else
        sum_filename = ['Sum_growth_',num2str(i),'.TXT'];
    end
    
    AEfilename = [Address_root,sum_filename];
    Sum_data = importdata(AEfilename,' ',4);
    
    SumData_temp = Sum_data.data;
    
    if i==1
        SumData = SumData_temp;
    else
        SumData(1,i+1) = SumData(1,i) + SumData_temp(1,2);
        SumData(2,i+1) = SumData(2,i) + SumData_temp(2,2);
    end
    
end

% Specify BinWidth for computing waveform entropy
BinWidth1 = 0.01;
BinWidth2 = 0.001;
BinWidth3 = 0.0001;

% Generate variables with hit filter data
Ch1dtf = HitData_Ch1.data;
Ch2dtf = HitData_Ch2.data;

Ch1kd = Ch1dtf(:,18);
Ch2kd = Ch2dtf(:,18);

% In this for loop, first AE entropy of load filter and Delta T filtered
% data are calculated. 
disp('AE entropy of waveforms for load and Delta T filtered hits are being calculated ...')

y1=ones(para_num,1);
y2=ones(para_num,1);
% For Channel 1 we have:
for k = 1:length(Ch1dtf)
    
    for i=1:para_num-1 % i=1 in test CTAA05, i=1:para_num-1 in rest
        
        if Ch1kd(k) < SumData(1,i+1)
            if i < 10
                Wave_filename_Ch1 = strcat('Growth - Waveform\Test_growth_0',num2str(i),'_1_',num2str(Ch1kd(k)),'.txt');
            else
                Wave_filename_Ch1 = strcat('Growth - Waveform\Test_growth_',num2str(i),'_1_',num2str(Ch1kd(k)),'.txt');
            end
            WaveAddress = [Address_root,Wave_filename_Ch1];
            Waveform_DTF_1_t = importdata(WaveAddress,' ',10);
        
            if Waveform_DTF_1_t.data - Ch1dtf(k,2) < -2e-5
                while Waveform_DTF_1_t.data - Ch1dtf(k,2) < -2e-5
                    if i < 10
                        Wave_filename_Ch1 = strcat('Growth - Waveform\Test_growth_0',num2str(i),'_1_',num2str(Ch1kd(k)+y1(i)),'.txt');
                    else
                        Wave_filename_Ch1 = strcat('Growth - Waveform\Test_growth_',num2str(i),'_1_',num2str(Ch1kd(k)+y1(i)),'.txt');
                    end
                    
                    WaveAddress = [Address_root,Wave_filename_Ch1];
                    
                    if Ch1kd(k)+y1(i) > SumData(1,i+1)
                        disp(['Hit number ',num2str(Ch1kd(k)),' in Channel 1 is skipped ! No waveform is found for it !'])
                        break
                    end
                    
                    Waveform_DTF_1_t = importdata(WaveAddress,' ',10);
                    y1(i)=y1(i)+1;
                    
                end
                
                y1(i)=y1(i)-1;

                if i < 10
                    disp(['Waveform ',Wave_filename_Ch1,' is used instead of',' Test_0',num2str(i),'_1_',num2str(Ch1kd(k)),'.'])
                else
                    disp(['Waveform ',Wave_filename_Ch1,' is used instead of',' Test_',num2str(i),'_1_',num2str(Ch1kd(k)),'.'])
                end
                
            end
        
            if Waveform_DTF_1_t.data - Ch1dtf(k,2) > 2e-5
                while Waveform_DTF_1_t.data - Ch1dtf(k,2) > 2e-5
                    if i < 10
                        Wave_filename_Ch1 = strcat('Growth - Waveform\Test_growth_0',num2str(i),'_1_',num2str(Ch1kd(k)-y2(i)),'.txt');
                    else
                        Wave_filename_Ch1 = strcat('Growth - Waveform\Test_growth_',num2str(i),'_1_',num2str(Ch1kd(k)-y2(i)),'.txt');
                    end
                    WaveAddress = [Address_root,Wave_filename_Ch1];
                    
                    if Ch1kd(k)-y2(i) <= 0
                        disp(['Hit number ',num2str(Ch1kd(k)),' in Channel 1 is skipped ! No waveform is found for it !'])
                        break                    
                    end
                    
                    Waveform_DTF_1_t = importdata(WaveAddress,' ',10);
                    y2(i)=y2(i)+1;
                    
                end
                
                y2(i)=y2(i)-1;
            
                if i < 10
                    disp(['Waveform ',Wave_filename_Ch1,' is used instead of',' Test_0',num2str(i),'_1_',num2str(Ch1kd(k)),'.'])
                else
                    disp(['Waveform ',Wave_filename_Ch1,' is used instead of',' Test_',num2str(i),'_1_',num2str(Ch1kd(k)),'.'])
                end
        
            end
        
            break % This break makes the code to jump into waveform calculation after finding right waveform file

        elseif SumData(1,i+1) <= Ch1kd(k) && Ch1kd(k) < SumData(1,i+2)
        
            if i < 9
                Wave_filename_Ch1 = strcat('Growth - Waveform\Test_growth_0',num2str(i+1),'_1_',num2str(Ch1kd(k)-SumData(1,i+1)),'.txt');
            else 
                Wave_filename_Ch1 = strcat('Growth - Waveform\Test_growth_',num2str(i+1),'_1_',num2str(Ch1kd(k)-SumData(1,i+1)),'.txt');
            end
            
            if SumData(1,i+1) == Ch1kd(k)
                if i < 9
                    Wave_filename_Ch1 = strcat('Growth - Waveform\Test_growth_0',num2str(i+1),'_1_1.txt');
                else 
                    Wave_filename_Ch1 = strcat('Growth - Waveform\Test_growth_',num2str(i+1),'_1_1.txt');
                end
            end
            
            WaveAddress = [Address_root,Wave_filename_Ch1];
            Waveform_DTF_1_t = importdata(WaveAddress,' ',10);
        
            if Waveform_DTF_1_t.data - Ch1dtf(k,2) < -2e-5
                while Waveform_DTF_1_t.data - Ch1dtf(k,2) < -2e-5
                    if i < 9
                        Wave_filename_Ch1 = strcat('Growth - Waveform\Test_growth_0',num2str(i+1),'_1_',num2str(Ch1kd(k)+y1(i+1)-SumData(1,i+1)),'.txt');
                    else
                        Wave_filename_Ch1 = strcat('Growth - Waveform\Test_growth_',num2str(i+1),'_1_',num2str(Ch1kd(k)+y1(i+1)-SumData(1,i+1)),'.txt');
                    end
                    WaveAddress = [Address_root,Wave_filename_Ch1];
                    
                    if Ch1kd(k)+y1(i+1) > SumData(1,i+2)
                        disp(['Hit number ',num2str(Ch1kd(k)),' in Channel 1 is skipped ! No waveform is found for it !'])
                        break
                    end
                    
                    Waveform_DTF_1_t = importdata(WaveAddress,' ',10);
                    y1(i+1)=y1(i+1)+1;
                    
                end
                
                y1(i+1)=y1(i+1)-1;
            
                if i < 9
                    disp(['Waveform ',Wave_filename_Ch1,' is used instead of',' Test_0',num2str(i+1),'_1_',num2str(Ch1kd(k)-SumData(1,i+1)),'.'])
                else
                    disp(['Waveform ',Wave_filename_Ch1,' is used instead of',' Test_',num2str(i+1),'_1_',num2str(Ch1kd(k)-SumData(1,i+1)),'.'])
                end

            end
        
            if Waveform_DTF_1_t.data - Ch1dtf(k,2) > 2e-5
                while Waveform_DTF_1_t.data - Ch1dtf(k,2) > 2e-5
                    if i < 9
                        Wave_filename_Ch1 = strcat('Growth - Waveform\Test_growth_0',num2str(i+1),'_1_',num2str(Ch1kd(k)-y2(i+1)-SumData(1,i+1)),'.txt');
                    else
                        Wave_filename_Ch1 = strcat('Growth - Waveform\Test_growth_',num2str(i+1),'_1_',num2str(Ch1kd(k)-y2(i+1)-SumData(1,i+1)),'.txt');
                    end
                    WaveAddress = [Address_root,Wave_filename_Ch1];
                    
                    if Ch1kd(k)-y2(i+1) <= 0
                        disp(['Hit number ',num2str(k),' in Channel 1 is skipped ! No waveform is found for it !'])
                        break                    
                    end
                    
                    Waveform_DTF_1_t = importdata(WaveAddress,' ',10);
                    y2(i+1)=y2(i+1)+1;
                    
                end
                
                y2(i+1)=y2(i+1)-1;
            
                if i < 9
                    disp(['Waveform ',Wave_filename_Ch1,' is used instead of',' Test_0',num2str(i+1),'_1_',num2str(Ch1kd(k)-SumData(1,i+1)),'.'])
                else
                    disp(['Waveform ',Wave_filename_Ch1,' is used instead of',' Test_',num2str(i+1),'_1_',num2str(Ch1kd(k)-SumData(1,i+1)),'.'])
                end

            end
               
            break % This break makes the code to jump into waveform calculation after finding right waveform file
        
        end

    end
    
    if abs(Waveform_DTF_1_t.data - Ch1dtf(k,2)) > 2e-5
        disp(['Hit number ',num2str(Ch1kd(k)),' in Ch1dtf does not have any waveform !'])
        continue
    end
    
%     waveform file is imported.
    Waveform_DTF_1 = importdata(WaveAddress,' ',12);
    
%    values regarding hit lockout time are all zero. Those values are
%    removed in this step. The criteria for removing is if average value
%    from a specific point to the end is zero, then take the index number
%    and cut the rest out.
       
     for i=1:length(Waveform_DTF_1.data)
         m = mean(Waveform_DTF_1.data(i:length(Waveform_DTF_1.data)));
         if m == 0
             HitLockOutIndex=i;
             break
         elseif i==10240
             HitLockOutIndex=i;
         end
     end
     
     if HitLockOutIndex ~= length(Waveform_DTF_1.data)
         Waveform_DTF_1.data(HitLockOutIndex:length(Waveform_DTF_1.data),:)=[];
     end
          
%     Probability distribution of of data with binwidth of BinWidth value is
%     generated
      [prob_1,edge_1] = histcounts(Waveform_DTF_1.data,'Binwidth',BinWidth1,'Normalization','probability');
      [prob_2,edge_2] = histcounts(Waveform_DTF_1.data,'Binwidth',BinWidth2,'Normalization','probability');
      [prob_3,edge_3] = histcounts(Waveform_DTF_1.data,'Binwidth',BinWidth3,'Normalization','probability');
    
%     The variable Entropy contains 3 columns, column 1 stores probability
%     value, column 2 log of the probability, and column 3 multipication of
%     column 1 and 2.
      entdtf1_1 = zeros(max(max(length(prob_1),3)),3);
      entdtf1_1(1:length(prob_1),1) = prob_1.';
      
      for i = 1:length(entdtf1_1)
          if entdtf1_1(i,1) == 0
              entdtf1_1(i,2) = 0;
          else
              entdtf1_1(i,2) = log(entdtf1_1(i,1));
          end
      end
      
      for i = 1:length(entdtf1_1)
          entdtf1_1(i,3) = -entdtf1_1(i,1).*entdtf1_1(i,2);
      end
    
%     AE entropy is calculated by sumation of all entries of column 3
    Ent_dtf_Ch1_1(k,1) = sum(entdtf1_1(:,3)); % Keep Entropy value in column 1
    Ent_dtf_Ch1_1(k,2) = Ch1dtf(k,19);      % Keep Cycle value in column 2
    Ent_dtf_Ch1_1(k,3) = Ch1dtf(k,20);      % Keep corrected time value in column 3
    
%     Repeat the same actions for BinWidth2
      entdtf1_2 = zeros(length(prob_2),3);
      entdtf1_2(:,1) = prob_2.';
      
      for i = 1:length(entdtf1_2)
          if entdtf1_2(i,1) == 0
              entdtf1_2(i,2) = 0;
          else
              entdtf1_2(i,2) = log(entdtf1_2(i,1));
          end
      end
      
      for i = 1:length(entdtf1_2)
          entdtf1_2(i,3) = -entdtf1_2(i,1).*entdtf1_2(i,2);
      end
    
%     AE entropy is calculated by sumation of all entries of column 3
    Ent_dtf_Ch1_2(k,1) = sum(entdtf1_2(:,3)); % Keep Entropy value in column 1
    Ent_dtf_Ch1_2(k,2) = Ch1dtf(k,19);      % Keep Cycle value in column 2
    Ent_dtf_Ch1_2(k,3) = Ch1dtf(k,20);      % Keep corrected time value in column 3
    
%     Repeat the same actions for BinWidth3
      entdtf1_3 = zeros(length(prob_3),3);
      entdtf1_3(:,1) = prob_3.';
      
      for i = 1:length(entdtf1_3)
          if entdtf1_3(i,1) == 0
              entdtf1_3(i,2) = 0;
          else
              entdtf1_3(i,2) = log(entdtf1_3(i,1));
          end
      end
      
      for i = 1:length(entdtf1_3)
          entdtf1_3(i,3) = -entdtf1_3(i,1).*entdtf1_3(i,2);
      end
    
%     AE entropy is calculated by sumation of all entries of column 3
    Ent_dtf_Ch1_3(k,1) = sum(entdtf1_3(:,3)); % Keep Entropy value in column 1
    Ent_dtf_Ch1_3(k,2) = Ch1dtf(k,19);      % Keep Cycle value in column 2
    Ent_dtf_Ch1_3(k,3) = Ch1dtf(k,20);      % Keep corrected time value in column 3

end

% Find unique cycle numbers 
Unique_Cycle_Ch1(:,1) = unique(Ent_dtf_Ch1_1(:,2));

% Find average entropy value for unique cycle numbers for BinWidth1
Mean_Ent_Ch1_1 = zeros(length(Unique_Cycle_Ch1),1);
for i=1:length(Unique_Cycle_Ch1)
        Mean_Ent_Ch1_1(i,:) = mean(Ent_dtf_Ch1_1(Ent_dtf_Ch1_1(:,2)==Unique_Cycle_Ch1(i),1));
end

% Add average entropy and unique cycles to a new variable
Ent_dtf_Ch1_1_mean = [Mean_Ent_Ch1_1 Unique_Cycle_Ch1];

% Find average entropy value for unique cycle numbers for BinWidth2
Mean_Ent_Ch1_2 = zeros(length(Unique_Cycle_Ch1),1);
for i=1:length(Unique_Cycle_Ch1)
        Mean_Ent_Ch1_2(i,:) = mean(Ent_dtf_Ch1_2(Ent_dtf_Ch1_2(:,2)==Unique_Cycle_Ch1(i),1));
end

% Add average entropy and unique cycles to a new variable
Ent_dtf_Ch1_2_mean = [Mean_Ent_Ch1_2 Unique_Cycle_Ch1];

% Find average entropy value for unique cycle numbers for BinWidth3
Mean_Ent_Ch1_3 = zeros(length(Unique_Cycle_Ch1),1);
for i=1:length(Unique_Cycle_Ch1)
        Mean_Ent_Ch1_3(i,:) = mean(Ent_dtf_Ch1_3(Ent_dtf_Ch1_3(:,2)==Unique_Cycle_Ch1(i),1));
end

% Add average entropy and unique cycles to a new variable
Ent_dtf_Ch1_3_mean = [Mean_Ent_Ch1_3 Unique_Cycle_Ch1];

% Draw AE entropy results for Channel 1 BinWidth2 for the whole test
f5 = figure('units','normalized','outerposition',[0 0 0.6 0.7]);
plot(Ent_dtf_Ch1_2(:,2),Ent_dtf_Ch1_2(:,1),'ro',...
    Ent_dtf_Ch1_2_mean(:,2),Ent_dtf_Ch1_2_mean(:,1),'b-')
title (['Ch1: load and Delta T filtered hits - BW=',num2str(BinWidth2)])
xlabel('Cycle')
ylabel('Entropy')
% axis(limit_ent)
legend('Entropy values','Change of Average','Location','northwest','Orientation','vertical')
print(f5,'06_Ent_Ch1_BW_2.tiff','-dtiffn','-r300')

% Compare BinWidth Change on the results from Ch1
f10 = figure('units','normalized','outerposition',[0 0 0.8 0.9]);
plot(Ent_dtf_Ch1_1_mean(:,2),Ent_dtf_Ch1_1_mean(:,1),'ro',...
    Ent_dtf_Ch1_1_mean(:,2),Ent_dtf_Ch1_1_mean(:,1),'r-',...
    Ent_dtf_Ch1_2_mean(:,2),Ent_dtf_Ch1_2_mean(:,1),'bo',...
    Ent_dtf_Ch1_2_mean(:,2),Ent_dtf_Ch1_2_mean(:,1),'b-',...
    Ent_dtf_Ch1_3_mean(:,2),Ent_dtf_Ch1_3_mean(:,1),'go',...
    Ent_dtf_Ch1_3_mean(:,2),Ent_dtf_Ch1_3_mean(:,1),'g-')
title ('Effect of BW Change on AE Entropy in Channel 1')
xlabel('Cycle')
ylabel('Entropy')
% axis(limit_ent)
legend('Average Entropy (BW 0.01)','Change of Average','Average Entropy (BW 0.001)',...
    'Change of Average','Average Entropy (BW 0.0001','Change of Average',...
    'Location','best','Orientation','vertical')
print(f10,'10_Ent_Ch1_BW_Compare.tiff','-dtiffn','-r300')

disp ('------------------------------------------------------------------')
disp('Channel 1 data is completed, starting Channel 2 ...')

t1=ones(para_num,1);
t2=ones(para_num,1);
for k = 1:length(Ch2dtf)
    
    for i=1:para_num-1 %1:para_num-1
        
        if Ch2kd(k) < SumData(2,i+1)
            if i < 10
                Wave_filename_Ch2 = strcat('Growth - Waveform\Test_growth_0',num2str(i),'_2_',num2str(Ch2kd(k)),'.txt');
            else
                Wave_filename_Ch2 = strcat('Growth - Waveform\Test_growth_',num2str(i),'_2_',num2str(Ch2kd(k)),'.txt');
            end
            WaveAddress = [Address_root,Wave_filename_Ch2];
            Waveform_DTF_2_t = importdata(WaveAddress,' ',10);
        
            if Waveform_DTF_2_t.data - Ch2dtf(k,2) < -2e-5
                while Waveform_DTF_2_t.data - Ch2dtf(k,2) < -2e-5
                    if i < 10
                        Wave_filename_Ch2 = strcat('Growth - Waveform\Test_growth_0',num2str(i),'_2_',num2str(Ch2kd(k)+t1(i)),'.txt');
                    else
                        Wave_filename_Ch2 = strcat('Growth - Waveform\Test_growth_',num2str(i),'_2_',num2str(Ch2kd(k)+t1(i)),'.txt');
                    end
                    WaveAddress = [Address_root,Wave_filename_Ch2];
                    
                    if Ch2kd(k)+t1(i) > SumData(2,i+1)
                        disp(['Hit number ',num2str(k),' in Channel 2 is skipped ! No waveform is found for it !'])
                        break
                    end
                    
                    Waveform_DTF_2_t = importdata(WaveAddress,' ',10);
                    t1(i)=t1(i)+1;
                    
                end
                
                t1(i)=t1(i)-1;

                if i < 10
                    disp(['Waveform ',Wave_filename_Ch2,' is used instead of',' Test_0',num2str(i),'_2_',num2str(Ch2kd(k)),'.'])
                else
                    disp(['Waveform ',Wave_filename_Ch2,' is used instead of',' Test_',num2str(i),'_2_',num2str(Ch2kd(k)),'.'])
                end
                
            end
        
            if Waveform_DTF_2_t.data - Ch2dtf(k,2) > 2e-5
                while Waveform_DTF_2_t.data - Ch2dtf(k,2) > 2e-5
                    if i < 10
                        Wave_filename_Ch2 = strcat('Growth - Waveform\Test_growth_0',num2str(i),'_2_',num2str(Ch2kd(k)-t2(i)),'.txt');
                    else
                        Wave_filename_Ch2 = strcat('Growth - Waveform\Test_growth_',num2str(i),'_2_',num2str(Ch2kd(k)-t2(i)),'.txt');
                    end
                    WaveAddress = [Address_root,Wave_filename_Ch2];
                    
                    if Ch2kd(k)-t2(i) <= 0
                        disp(['Hit number ',num2str(Ch2kd(k)),' in Channel 2 is skipped ! No waveform is found for it !'])
                        break                    
                    end
                    
                    Waveform_DTF_2_t = importdata(WaveAddress,' ',10);
                    t2(i)=t2(i)+1;
                    
                end
                
                t2(i)=t2(i)-1;
            
                if i < 10
                    disp(['Waveform ',Wave_filename_Ch2,' is used instead of',' Test_0',num2str(i),'_2_',num2str(Ch2kd(k)),'.'])
                else
                    disp(['Waveform ',Wave_filename_Ch2,' is used instead of',' Test_',num2str(i),'_2_',num2str(Ch2kd(k)),'.'])
                end
        
            end
        
            break % This break makes the code to jump into waveform calculation after finding right waveform file

        elseif SumData(2,i+1) <= Ch2kd(k) && Ch2kd(k) < SumData(2,i+2)
        
            if i < 9
                Wave_filename_Ch2 = strcat('Growth - Waveform\Test_growth_0',num2str(i+1),'_2_',num2str(Ch2kd(k)-SumData(2,i+1)),'.txt');
            else 
                Wave_filename_Ch2 = strcat('Growth - Waveform\Test_growth_',num2str(i+1),'_2_',num2str(Ch2kd(k)-SumData(2,i+1)),'.txt');
            end
            
            if SumData(2,i+1) == Ch2kd(k)
                if i < 9
                    Wave_filename_Ch2 = strcat('Growth - Waveform\Test_growth_0',num2str(i+1),'_2_1.txt');
                else 
                    Wave_filename_Ch2 = strcat('Growth - Waveform\Test_growth_',num2str(i+1),'_2_1.txt');
                end
            end
            
            WaveAddress = [Address_root,Wave_filename_Ch2];
            Waveform_DTF_2_t = importdata(WaveAddress,' ',10);
        
            if Waveform_DTF_2_t.data - Ch2dtf(k,2) < -2e-5
                while Waveform_DTF_2_t.data - Ch2dtf(k,2) < -2e-5
                    if i < 9
                        Wave_filename_Ch2 = strcat('Growth - Waveform\Test_growth_0',num2str(i+1),'_2_',num2str(Ch2kd(k)+t1(i+1)-SumData(2,i+1)),'.txt');
                    else
                        Wave_filename_Ch2 = strcat('Growth - Waveform\Test_growth_',num2str(i+1),'_2_',num2str(Ch2kd(k)+t1(i+1)-SumData(2,i+1)),'.txt');
                    end
                    
                    WaveAddress = [Address_root,Wave_filename_Ch2];
                    
                    if Ch2kd(k)+t1(i+1) > SumData(2,i+2)
                        disp(['Hit number ',num2str(Ch2kd(k)),' in Channel 2 is skipped ! No waveform is found for it !'])
                        break
                    end
                    
                    Waveform_DTF_2_t = importdata(WaveAddress,' ',10);
                    t1(i+1)=t1(i+1)+1;
                    
                end
                
                t1(i+1)=t1(i+1)-1;
            
                if i < 9
                    disp(['Waveform ',Wave_filename_Ch2,' is used instead of',' Test_0',num2str(i+1),'_2_',num2str(Ch2kd(k)-SumData(2,i+1)),'.'])
                else
                    disp(['Waveform ',Wave_filename_Ch2,' is used instead of',' Test_',num2str(i+1),'_2_',num2str(Ch2kd(k)-SumData(2,i+1)),'.'])
                end

            end
        
            if Waveform_DTF_2_t.data - Ch2dtf(k,2) > 2e-5
                while Waveform_DTF_2_t.data - Ch2dtf(k,2) > 2e-5
                    if i < 9
                        Wave_filename_Ch2 = strcat('Growth - Waveform\Test_growth_0',num2str(i+1),'_2_',num2str(Ch2kd(k)-SumData(2,i+1)-t2(i+1)),'.txt');
                    else
                        Wave_filename_Ch2 = strcat('Growth - Waveform\Test_growth_',num2str(i+1),'_2_',num2str(Ch2kd(k)-SumData(2,i+1)-t2(i+1)),'.txt');
                    end
                    WaveAddress = [Address_root,Wave_filename_Ch2];
                   
                    if Ch2kd(k)-t2(i+1) <= 0
                        disp(['Hit number ',num2str(Ch2kd(k)),' in Channel 2 is skipped ! No waveform is found for it !'])
                        break                    
                    end
                    
                    Waveform_DTF_2_t = importdata(WaveAddress,' ',10);
                    t2(i+1)=t2(i+1)+1;
                    
                end
                
                t2(i+1)=t2(i+1)-1;
            
                if i < 9
                    disp(['Waveform ',Wave_filename_Ch2,' is used instead of',' Test_0',num2str(i+1),'_2_',num2str(Ch2kd(k)-SumData(2,i+1)),'.'])
                else
                    disp(['Waveform ',Wave_filename_Ch2,' is used instead of',' Test_',num2str(i+1),'_2_',num2str(Ch2kd(k)-SumData(2,i+1)),'.'])
                end

            end
        
            break % This break makes the code to jump into waveform calculation after finding right waveform file
        
        end
       
    end 
    
    if abs(Waveform_DTF_2_t.data - Ch2dtf(k,2)) > 2e-5
       disp(['Hit number ',num2str(Ch2kd(k)),' in Ch2dtf does not have any waveform !'])
       continue
    end
        
    
%     waveform file is imported.
    Waveform_DTF_2 = importdata(WaveAddress,' ',12);
    
%    values regarding hit lockout time are all zero. Those values are
%    removed in this step. The criteria for removing is if average value
%    from a specific point to the end is zero, then take the index number
%    and cut the rest out.
       
     for i=1:length(Waveform_DTF_2.data)
         m = mean(Waveform_DTF_2.data(i:length(Waveform_DTF_2.data)));
         if m == 0
             HitLockOutIndex=i;
             break
         elseif i==10240
             HitLockOutIndex=i;
         end
      end
     
      if HitLockOutIndex ~= length(Waveform_DTF_2.data)
         Waveform_DTF_2.data(HitLockOutIndex:length(Waveform_DTF_2.data),:)=[];
      end
          
%     Probability distribution of of data with binwidth of BinWidth value is
%     generated
      [prob_1,edge_1] = histcounts(Waveform_DTF_2.data,'Binwidth',BinWidth1,'Normalization','probability');
      [prob_2,edge_2] = histcounts(Waveform_DTF_2.data,'Binwidth',BinWidth2,'Normalization','probability');
      [prob_3,edge_3] = histcounts(Waveform_DTF_2.data,'Binwidth',BinWidth3,'Normalization','probability');
    
%     The variable Entropy contains 3 columns, column 1 stores probability
%     value, column 2 log of the probability, and column 3 multipication of
%     column 1 and 2.
      entdtf2_1 = zeros(max(max(length(prob_1),3)),3);
      entdtf2_1(1:length(prob_1),1) = prob_1.';
      
      for i = 1:length(entdtf2_1)
          if entdtf2_1(i,1) == 0
              entdtf2_1(i,2) = 0;
          else
              entdtf2_1(i,2) = log(entdtf2_1(i,1));
          end
      end
      
      for i = 1:length(entdtf2_1)
          entdtf2_1(i,3) = -entdtf2_1(i,1).*entdtf2_1(i,2);
      end
    
%     AE entropy is calculated by sumation of all entries of column 3
      Ent_dtf_Ch2_1(k,1) = sum(entdtf2_1(:,3)); % Keep Entropy value in column 1
      Ent_dtf_Ch2_1(k,2) = Ch2dtf(k,19);      % Keep Cycle value in column 2
      Ent_dtf_Ch2_1(k,3) = Ch2dtf(k,20);      % Keep corrected time value in column 3
    
%     Repeat the same actions for BinWidth2
      entdtf2_2 = zeros(length(prob_2),3);
      entdtf2_2(:,1) = prob_2.';
      
      for i = 1:length(entdtf2_2)
          if entdtf2_2(i,1) == 0
              entdtf2_2(i,2) = 0;
          else
              entdtf2_2(i,2) = log(entdtf2_2(i,1));
          end
      end
      
      for i = 1:length(entdtf2_2)
          entdtf2_2(i,3) = -entdtf2_2(i,1).*entdtf2_2(i,2);
      end
    
%     AE entropy is calculated by sumation of all entries of column 3
      Ent_dtf_Ch2_2(k,1) = sum(entdtf2_2(:,3)); % Keep Entropy value in column 1
      Ent_dtf_Ch2_2(k,2) = Ch2dtf(k,19);      % Keep Cycle value in column 2
      Ent_dtf_Ch2_2(k,3) = Ch2dtf(k,20);      % Keep corrected time value in column 3
    
%     Repeat the same actions for BinWidth3
      entdtf2_3 = zeros(length(prob_3),3);
      entdtf2_3(:,1) = prob_3.';
      
      for i = 1:length(entdtf2_3)
          if entdtf2_3(i,1) == 0
              entdtf2_3(i,2) = 0;
          else
              entdtf2_3(i,2) = log(entdtf2_3(i,1));
          end
      end
      
      for i = 1:length(entdtf2_3)
          entdtf2_3(i,3) = -entdtf2_3(i,1).*entdtf2_3(i,2);
      end
    
%     AE entropy is calculated by sumation of all entries of column 3
     Ent_dtf_Ch2_3(k,1) = sum(entdtf2_3(:,3)); % Keep Entropy value in column 1
     Ent_dtf_Ch2_3(k,2) = Ch2dtf(k,19);      % Keep Cycle value in column 2
     Ent_dtf_Ch2_3(k,3) = Ch2dtf(k,20);      % Keep corrected time value in column 3

end

% Find unique cycle numbers 
Unique_Cycle_Ch2(:,1) = unique(Ent_dtf_Ch2_1(:,2));

% Find average entropy value for unique cycle numbers for BinWidth2
Mean_Ent_Ch2_1 = zeros(length(Unique_Cycle_Ch2),1);
for i=1:length(Unique_Cycle_Ch2)
        Mean_Ent_Ch2_1(i,:) = mean(Ent_dtf_Ch2_1(Ent_dtf_Ch2_1(:,2)==Unique_Cycle_Ch2(i),1));
end

% Add average entropy and unique cycles to a new variable
Ent_dtf_Ch2_1_mean = [Mean_Ent_Ch2_1 Unique_Cycle_Ch2];

% Find average entropy value for unique cycle numbers for BinWidth1
Mean_Ent_Ch2_2 = zeros(length(Unique_Cycle_Ch2),1);
for i=1:length(Unique_Cycle_Ch2)
        Mean_Ent_Ch2_2(i,:) = mean(Ent_dtf_Ch2_2(Ent_dtf_Ch2_2(:,2)==Unique_Cycle_Ch2(i),1));
end

% Add average entropy and unique cycles to a new variable
Ent_dtf_Ch2_2_mean = [Mean_Ent_Ch2_2 Unique_Cycle_Ch2];

% Find average entropy value for unique cycle numbers for BinWidth1
Mean_Ent_Ch2_3 = zeros(length(Unique_Cycle_Ch2),1);
for i=1:length(Unique_Cycle_Ch2)
        Mean_Ent_Ch2_3(i,:) = mean(Ent_dtf_Ch2_3(Ent_dtf_Ch2_3(:,2)==Unique_Cycle_Ch2(i),1));
end

% Add average entropy and unique cycles to a new variable
Ent_dtf_Ch2_3_mean = [Mean_Ent_Ch2_3 Unique_Cycle_Ch2];

disp('AE waveform entropy calculation of load and Delta T filtered data is completed.')
disp ('------------------------------------------------------------------')

% Draw AE entropy results for Channel 2 BinWidth2
f12 = figure('units','normalized','outerposition',[0 0 0.6 0.7]);
plot(Ent_dtf_Ch2_2(:,2),Ent_dtf_Ch2_2(:,1),'ro',...
    Ent_dtf_Ch2_2_mean(:,2),Ent_dtf_Ch2_2_mean(:,1),'b-')
title (['Ch2: load and Delta T filtered hits - BW=',num2str(BinWidth2)])
xlabel('Cycle')
ylabel('Entropy')
% axis(limit_ent)
legend('Entropy values','Change of Average','Location','northwest','Orientation','vertical')
print(f12,'13_Ent_Ch2_BW_2.tiff','-dtiffn','-r300')

% Compare BinWidth Change on the results from Ch2
f17 = figure('units','normalized','outerposition',[0 0 0.8 0.9]);
plot(Ent_dtf_Ch2_1_mean(:,2),Ent_dtf_Ch2_1_mean(:,1),'ro',...
    Ent_dtf_Ch2_1_mean(:,2),Ent_dtf_Ch2_1_mean(:,1),'r-',...
    Ent_dtf_Ch2_2_mean(:,2),Ent_dtf_Ch2_2_mean(:,1),'bo',...
    Ent_dtf_Ch2_2_mean(:,2),Ent_dtf_Ch2_2_mean(:,1),'b-',...
    Ent_dtf_Ch2_3_mean(:,2),Ent_dtf_Ch2_3_mean(:,1),'go',...
    Ent_dtf_Ch2_3_mean(:,2),Ent_dtf_Ch2_3_mean(:,1),'g-')
title ('Effect of BW Change on AE Entropy in Channel 2')
xlabel('Cycle')
ylabel('Entropy')
% axis(limit_ent)
legend('Average Entropy (BW 0.01)','Change of Average','Average Entropy (BW 0.001)',...
    'Change of Average','Average Entropy (BW 0.0001)','Change of Average',...
    'Location','best','Orientation','vertical')
print(f17,'17_Ent_Ch2_BW_Compare.tiff','-dtiffn','-r300')

% Compare mean entropy Change in Channel 1 and 2 for BinWidth2 (0.001)
f18 = figure('units','normalized','outerposition',[0 0 0.6 0.7]);
plot(Ent_dtf_Ch1_2_mean(:,2),Ent_dtf_Ch1_2_mean(:,1),'ro',...
    Ent_dtf_Ch1_2_mean(:,2),Ent_dtf_Ch1_2_mean(:,1),'r-',...
    Ent_dtf_Ch2_2_mean(:,2),Ent_dtf_Ch2_2_mean(:,1),'bx',...
    Ent_dtf_Ch2_2_mean(:,2),Ent_dtf_Ch2_2_mean(:,1),'b--')
title (['Comparison of AE Entropy in Ch 1 and 2, BW= ',num2str(BinWidth2)])
xlabel('Cycle')
ylabel('Entropy')
% axis(limit_ent)
legend('Average Entropy - Ch1','Change of Average','Average Entropy - Ch2',...
    'Change of Average','Location','southwest','Orientation','vertical')
print(f18,'18_Ent_Ch1_Ch2_BW_2_Compare.tiff','-dtiffn','-r300')

% % % % % % % Write entropy values of Ent_dtf_Ch1_2 and Ent_dtf_Ch2_1 to an excel file
filename1 = 'Ch1_BW_0.001_Ent_Cycle_Time.xlsx';
filename2 = 'Ch2_BW_0.001_Ent_Cycle_Time.xlsx';
xlswrite(filename1,Ent_dtf_Ch1_2)
xlswrite(filename2,Ent_dtf_Ch2_2)

t_total = toc;

disp('Summary:')
disp(['Test completed after ',num2str(max(max(Ch1dtf(:,19),Ch2dtf(:,19)))),' cycles.'])
disp(['Channel 1 RAW data contains ',num2str(SumData(1,1+para_num)),' hits.'])
disp(['Channel 2 RAW data contains ',num2str(SumData(2,1+para_num)),' hits.'])
disp([num2str(length(Ch2dtf)),' hits (from each Channel) satisfy load and Delta T conditions and used to calculate waveform entropy.'])
disp('Hits arrival time in Hit data file (Test_growth_##) is matChed with arrival time in waveform files.')
disp(['Bin widths used to calculate AE hit waveform entropies are: ',num2str(BinWidth1,'%.5g'),', ',num2str(BinWidth2,'%.5g'),', ',num2str(BinWidth3,'%.5g'),'.'])
disp('Waveform entropy data generated using bin width 0.001 for each Channel are exported.')
disp(['Total run time is ',num2str(t_total,'%.1f'),' seconds.'])
