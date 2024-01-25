%%% New formatting script
% data_folder = 'C:\Users\Somlab\Box\BensmaiaLab\ProjectFolders\DARPA\Data\RawData\Pinot\Electrode_31and41\ca_threshold';
data_folder = 'B:\ProjectFolders\DARPA\Data\RawData\Whistlepig\Electrodde_3and15\ThresholdTask\Cathodic';
file_list = dir(data_folder);


%% Load and parse .rsp files

ii = 1;
block_struct = struct();
for i = 1:length(file_list)
    % Check if file is .rsp
    if ~contains(file_list(i).name, '.rsp') || ~contains(file_list(i).name, 'Electrode')
        continue
    end
    
    % Parse filename
    fname_split = strsplit(file_list(i).name, '_');
    %data_struct(ii).StudyID = fname_split{1};
    block_struct(ii).Animal = fname_split{2}(1:end-7);
   % block_struct(ii).Protocol = fname_split{3};
    t_idx = find(fname_split{4} == 'T');
    block_struct(ii).Date = datestr(datenum(fname_split{4}(1:t_idx-1), 'yyyymmdd'), 1);
    %data_struct(ii).Time = datestr(datenum(fname_split{4}(t_idx+1:end-4), 'hhmmss'), 13);
    
     % Load the data
    temp_data = readcell(fullfile(data_folder, file_list(i).name), 'FileType', 'text', 'NumHeaderLines', 1);
    if size(temp_data,2) == 10
          % Remove 2nd column
         temp_data = temp_data(:,[1:3,5,6,8,10]);
    elseif size(temp_data,2) == 11
          temp_data = temp_data(:,[1:3,5,6,8,10]);
    end
%     
%     % Remove aborted trials
     abort_idx = strcmpi(temp_data(:,7), 'empty response') | strcmpi(temp_data(:,7), 'empty');
     temp_data = temp_data(~abort_idx,:); 
% %     % Convert to table
    response_table = cell2table(temp_data(1:end,:), 'VariableNames',{'Trial','CorrectAnswer', 'CorrectAnswerText', 'TestStimAmp', 'TestStimFreq','TestElectrodeNum','Response'});
    block_struct(ii).ResponseTable = response_table;



    ii = ii + 1;
end
%%


for b = 1:length(block_struct)
    num_tnum = max(block_struct(b).ResponseTable.Trial);
    correct_first = false(num_tnum,1);
    num_attempts = zeros(num_tnum,1);
    for t = 1:num_tnum
       idx = find(block_struct(b).ResponseTable.Trial == t, 1, 'first');
        if (strcmp(block_struct(b).ResponseTable.Response{idx}, 'correct'))
            correct_first(t) = true;
        end
        num_attempts(t) = length(idx);
    end   
    block_struct(b).percent_correct = sum(correct_first) / num_tnum;
    block_struct(b).percent_correct = block_struct(b).percent_correct * 100;
    block_struct(b).trials = num_tnum;
    
    block_struct(b).correct_trials = sum(correct_first);
    
%    plot(block_struct(b).percent_correct) 
end
%%
figure; hold on
x = [1: [length(block_struct)]];
plot(x, [block_struct.percent_correct], 'LineWidth', 3)
date_strs = '';
for t = 1:length(block_struct)

 date_strs{t} = block_struct(t).Date;

end
  ax = gca;
      ax.FontSize = 20;

%set(gca, 'XTick', [0:20], 'XTickLabels', date_strs')
xlabel('Session', 'FontSize', 18)
ylabel('Percent Correct', 'FontSize', 18)
%title('Whistlepig Electrode 23')

%% get pdetect and dprime 

for i = 1:length(block_struct)
    [u_test_amps, ~, ia] = unique(block_struct(i).ResponseTable.TestStimAmp);
    p_detect_comb = zeros([length(u_test_amps),1]);
    for j = 1:length(u_test_amps)
        correct_idx = strcmp(block_struct(i).ResponseTable.Response(ia == j), 'correct');
        p_detect_comb(j) = sum(correct_idx) / length(correct_idx);
        
    end
    % Correct for catch trials
    if any(u_test_amps == 0)
        catch_idx = find(u_test_amps == 0);
        p_detect_comb(catch_idx) = 1 - p_detect_comb(catch_idx);
    end
    
    % Compute d' from p(detect) where  d' = norminv(p(hit)) - norminv(p(false alarm))
    dprime = NaN([length(u_test_amps),1]);
    pmiss = p_detect_comb(1);
    if pmiss == 0 % Correct for 0 false alarm
        pmiss = 0.001;
    end
    for j = 1:length(dprime)-1
        phit = p_detect_comb(j+1);
        if phit == 1 % Correct for infinite hit rate
            phit = .999;
        end
        dprime(j+1) = norminv(phit) - norminv(pmiss);        
    end
    
    % Make a table & add to struct
    block_struct(i).DetectionRates = array2table([u_test_amps, p_detect_comb, dprime], 'VariableNames', {'Amplitude', 'pDetect', 'dPrime'});
end
%% plotting day by day
% line_color = winter(length(pre_idx))
c = ColorGradient(rgb(255, 235, 238), rgb(183, 28, 28), length(block_struct));
for i = 1:length(block_struct)
  % figure('Name', sprintf('%s - %s', block_struct(i).Animal, block_struct(i).Date));
    subplot(1,2,1); hold on

    plot(block_struct(i).DetectionRates{:,1}, block_struct(i).DetectionRates{:,2},'o-','MarkerSize', 5, 'Color', c(i,:), 'LineWidth', 3)
    
    ax = gca;
    ax.FontSize = 18;
    xlabel(sprintf('Amplitude (%sA)', GetUnicodeChar('mu')),'FontSize', 18 )
    ylabel('p(Detected)','FontSize',18)

    subplot(1,2,2); hold on
    
     plot(block_struct(i).DetectionRates{:,1}, block_struct(i).DetectionRates{:,3},'o-','MarkerSize', 5, 'Color', c(i,:), 'LineWidth', 3)
      ax = gca;
      ax.FontSize = 18;
   yline(1.35,'-', 'Threshold ', 'FontSize',18,'LabelHorizontalAlignment','right', 'LineWidth',2);
    
    xlabel(sprintf('Amplitude (%sA)', GetUnicodeChar('mu')),'FontSize', 18)
    ylabel('d''','FontSize',18) 
    box off
   
end


%% plotting all combined

comb_idx = [1:[length(block_struct)]];

bigtable = cat(1, block_struct(comb_idx).ResponseTable);

for i = 1:size(bigtable)
 [u_test_amps_comb, ~, ia] = unique(bigtable.TestStimAmp);
    p_detect_comb = zeros([length(u_test_amps_comb),1]);
    for j = 1:length(u_test_amps_comb)
        correct_idx_comb = strcmp(bigtable.Response(ia == j), 'correct');
        p_detect_comb(j) = sum(correct_idx_comb) / length(correct_idx_comb);
        
    end
    % Correct for catch trials
    if any(u_test_amps_comb == 0)
        catch_idx_comb = find(u_test_amps_comb == 0);
        p_detect_comb(catch_idx_comb) = 1 - p_detect_comb(catch_idx_comb);
    end
    
    % Compute d' from p(detect) where  d' = norminv(p(hit)) - norminv(p(false alarm))
    dprime_comb = NaN([length(u_test_amps_comb),1]);
    pmiss_comb = p_detect_comb(1);
    if pmiss_comb == 0 % Correct for 0 false alarm
        pmiss_comb = 0.001;
    end
    for j = 1:length(dprime_comb)-1
        phit_comb = p_detect_comb(j+1);
        if phit_comb == 1 % Correct for infinite hit rate
            phit_comb = .999;
        end
        dprime_comb(j+1) = norminv(phit_comb) - norminv(pmiss_comb);        
    end
    
    % Make a table & add to struct
    combtable_DetectionRates = array2table([u_test_amps_comb, p_detect_comb, dprime_comb], 'VariableNames', {'Amplitude', 'pDetect', 'dPrime'});   


end




%% 
 a = combtable_DetectionRates{:,3};

    a(isnan(a)) = 0;

%% plotting for bigtable
subplot(1,2,1);
hold on;
plot(combtable_DetectionRates{:,1}, combtable_DetectionRates{:,2},'o-', 'MarkerSize', 5,'Color', rgb(33, 33, 33), 'LineWidth', 4);
ax = gca;
ax.FontSize = 18;
xlabel(sprintf('Stimulus Amplitude (%sA)', GetUnicodeChar('mu')), 'FontSize', 18);
ylabel('p(Detected)', 'FontSize', 18);
axis square



subplot(1,2,2);
hold on;
plot(combtable_DetectionRates{:,1}, combtable_DetectionRates{:,3},'o-','MarkerSize', 5, 'Color', rgb(33, 33, 33), 'LineWidth', 4);
ax = gca;
ax.FontSize = 18;
xlabel(sprintf('Stimulus Amplitude (%sA)', GetUnicodeChar('mu')), 'FontSize', 18);
ylabel('d''', 'FontSize', 18);
y_line = 1.35; % Specify the y-value for the y-line
axis square
% Find the intersection

x_data = combtable_DetectionRates{:,1}; 
y_data = a;
x_intersection = interp1(y_data, x_data, y_line);
 fix_x_intersection = fix(x_intersection);
% yline(1.35,'-', sprintf('Threshold 24', GetUnicodeChar('mu')),'LabelHorizontalAlignment','left', 'FontSize',18);
  yline(1.35,'-', sprintf('Threshold %.0f %sA', fix_x_intersection, GetUnicodeChar('mu')),'LabelHorizontalAlignment','left', 'FontSize',18);
% Plot the intersection point
if ~isnan(x_intersection)
    plot(x_intersection, y_line);
%      text(x_intersection, y_line, sprintf('Threshold %.2f', x_intersection), 'FontSize', 18, 'VerticalAlignment', 'bottom');
    
else
    disp('No intersection found.');
end

%% charles new function for residuals and resnorm

y_data1 = a;
% y_data1 = combtable_DetectionRates.dPrime;
x_data1 = combtable_DetectionRates.Amplitude;
sigfun = @(c,x) (c(3) .* (1./(1 + exp(-c(1).*(x-c(2)))))) + c(4); 
x0 = [1.7/mean(a),...% c(1)
      mean(a),...% c(2)
      max(combtable_DetectionRates.pDetect),...% c(3)
      min(combtable_DetectionRates.pDetect)];... % c(4)

% [SigmoidFun, coeffs, rnorm, residuals, jnd, warning] = FitSigmoid(xdata,ydata, 'PlotFit', true, 'CoeffInit', [1,15,NaN,NaN], 'NumCoeffs', 3);
 coeffs = SearchSigmoid(x_data1, y_data1, [1, 18], true); 
[~, coeffs, rnorm, residuals, jnd, ~] = FitSigmoid(y_data1,x_data1, 'PlotFit', true, 'CoeffInit', [1,15,NaN,NaN], 'NumCoeffs', 3);

% [x, resnorm,residual] = lsqcurvefit(sigfun,x0, xdata,ydata,lb,ub);



%%

% option #1
wind_size = 50;
for i = 1:(size(bigtable,1) - wind_size + 1) % might be off by 1 trial
    current_window = i:i+wind_size-1;
    [detection_table{i}, coeff_table{i}] = AnalyzeDetectionTable(bigtable(current_window, :));
end

% Use bigtable or detection_table to figure out what the unique stim
% amplitudes are
% unique_stim_amp = ....; DONE
% detection_table_dprime = NaN(length(unique_stim_amp), DONE
% length(detection_table));
num_windows = size(detection_table, 2);
StimAmps = unique(bigtable.TestStimAmp); % actual stim amp values
num_stim_amp = length(StimAmps);
detection_table_dprimed = NaN(num_stim_amp, num_windows);



detection_table_dprimed_2 = [];
for j = 1:length(detection_table)
    % For each sliding window table, find what stim amps are there, and
    % then put the dprimes into the corresponding slot in the
    % detection_table_dprime
    for s = 1:size(detection_table{j}, 1)
        
        detection_table_dprimed_2(s,j) = detection_table{j}.dPrime(s);
    end
end



 detection_table_table = table(StimAmps, detection_table_dprimed_2);
for c = 1:size(detection_table_dprimed_2,2)
 xdata = detection_table_table.StimAmps(:,1);
 ydata = detection_table_dprimed_2(:,c);



[~, coeffs_more{c}, rnorm_more{c}, residuals_more{c}, jnd_more{c},~] = FitSigmoid(xdata, ydata);


end




hold on
quick = 1:length(rnorm_more);


for i = 1:length(rnorm_more)
    plot(i, rnorm_more{i}, '-o');
end

xlabel('Trials', 'FontSize', 18)
ylabel('Rsnorm', 'FontSize', 18)





%breaks computer; probably just continously looping with no end?
% hold on
% quick = 1:length(rnorm_more);
% 
% 
% for i = 1:length(rnorm_more)
%     plot(quick, rnorm_more{i}, 'o-');
% end


%wrong? no wrong!

% for t = 1:size(rnorm_more,2)
%     hold on
% 
%      quick = 1:10; %size(rnorm_more,2);
%      plot(quick, rnorm_more{quick},'o-')
% %     plot(quick, rnorm_more{quick}, 'o-')
% % xlabel = ('Trials');
% % ylabel = ('Rnorm');
% 
% end
% ran_try = randn(1001,1);
%  plot(quick, ran_try)

% [~, coeffs_more{c}, rnorm_more{c}, residuals_more{c}, jnd_more{c},~] = FitSigmoid(ydata, xdata);


 % [detection_table{i}, coeff_table{i}] = AnalyzeDetectionTable(bigtable(current_window, :));


% Now in the loop -pre alexandria - broken
% for j = 1:length(detection_table)
%     % For each sliding window table, find what stim amps are there, and
%     % then put the dprimes into the corresponding slot in the
%     % detection_table_dprime
% 
%     stim_amp_idx = find(StimAmps == detection_table{1,j}.StimAmp);
% 
%     for s = 1:size(detection_table{j}, 1)
% 
%         detection_table_dprimed(stim_amp_idx, j) = detection_table{j}.dPrime(s);
% 
%     end
% 
% 
% end


% % option 2
% wind_size = 20;
% for i = 1:(size(bigtable,1) - wind_size + 1) % might be off by 1 trial
%     current_window = i:i+wind_size-1;
%     [det_tmp, coeff_tmp] = AnalyzeDetectionTable(bigtable(current_window, :));
% 
% detection_table_dprime(:, i) = det_tmp.dPrime(:); %detection_table{i}.dPrime(:);
% end
% %     
% for i = 1:size(bigtable, 1)
% % for i = 1:length(block_struct)
% %     num_trials = size(bigtable,1);
%      num_tnum = max(block_struct(i).ResponseTable.Trial);
%         block_struct(i).trials = num_tnum;
%     num_trials = block_struct(i).trials;
%     max_idx = num_trials - wind_size;
% % %     running_perf = zeros(max_idx,1);
% %     detection_table= zeros(max_idx,1);
% %     coeff_table = zeros(max_idx,1);
%      for j = 1:max_idx
% 
% %            current_window = bigtable(j:j+wind_size-1);
% 
% %      [detection_table{j}, coeff_table{j}] = AnalyzeDetectionTable((block_struct().ResponseTable(:,:),(j:j+wind_size-1));
%        
% 
     % detection_table_pdetect(stim_amp_idx, j) = detection_table{j}.pDetect(s);

     %% sliding window for dprime
     %increase to 100 and moving through increaments of 10

     wind_size = 100;
for i = 1:(size(bigtable,1) - wind_size + 1) % might be off by 1 trial
    current_window = i:i+wind_size - 1;
    % [detection_table{i}, coeff_table{i}] = AnalyzeDetectionTable(bigtable(current_window, :));
end

%      wind_size = 100;
% for i = 1:(size(bigtable,1) - wind_size + 1) % might be off by 1 trial
%      current_window = i:i+wind_size-1;
%     next = i+10:i+wind_size-1;
%      [detection_table_next{i}, coeff_table_next{i}] = AnalyzeDetectionTable(bigtable(next, :));
%      [detection_table_current{i}, coeff_table_current{i}] = AnalyzeDetectionTable(bigtable(current_window, :));
% end



% temp = 1:1000;
% i = 10:10:1000

% i = 10:10:1000;
% for i = 10:10:100
% length(temp(end-i:end))
% end
% 
% i = 10:10:100;
% for r = 1:length(i)
% length(temp(end-i(r):end))
% end

% inc = 20;
% i = 10:inc:100
% inc = 50;

% i = 10:inc:100


 %% starting over plotting
% 
% dprime_threshold = 1.35;
% % sigfun = @(c,x) (1./(1 + exp(-c(1).*(x-c(2)))));
% 
% 
% for i = 1:length(block_struct)
% 
%     x = block_struct(i).DetectionRates.Amplitude;
%     y = block_struct(i).DetectionRates.dPrime;
%     scatter(x,y,50, rgb(33, 33, 33), "filled")
%     plot(x,y, 'Color', rgb(33, 33, 33), 'LineStyle', ':')
% 
%     xq = linspace(x(1), x(end));
%     yq = sigfun(block_struct(i).DetectionRates.dPrime, xq);
%     f= fit(x,y);
% %     plot(xq, yq, 'Color', rgb(33, 33, 33))
% %     [~,b] = min(abs(yq-dprime_threshold));
% %     plot([0 xq(b) xq(b)], [dprime_threshold, dprime_threshold -1], 'Color', [.4 .4 .4], 'LineStyle','--')
% end
% 
% 


% box off
% title("Electrode 12")
% 'o-', 'MarkerFaceColor'
% %% trying to find poin of threshold 
% 
% dprime_threshold = 1.35;
% sigfun = @(x,c) (c(1) .* (1./(1 + exp(-c(1).*(x-c(2)))))) + c(2);
% 
% for i = 1:length(block_struct)
% 
%     subplot(1,2,1); hold on
% 
%     plot(block_struct(i).DetectionRates{:,1}, block_struct(i).DetectionRates{:,2},'o-',  'LineWidth', 3)
% 
%     ax = gca;
%     ax.FontSize = 18;
%     xlabel(sprintf('Amplitude (%sA)', GetUnicodeChar('mu')),'FontSize', 18 )
%     ylabel('p(Detected)','FontSize',18)
% 
%     subplot(1,2,2); hold on
%     
%      plot(block_struct(i).DetectionRates{:,1}, block_struct(i).DetectionRates{:,3},'o-', 'LineWidth', 3)
%       ax = gca;
%       ax.FontSize = 18;
%    yline(1.35,'-', 'Threshold ', 'FontSize',18,'LabelHorizontalAlignment','left', 'LineWidth',2);
%     
%     xlabel(sprintf('Amplitude (%sA)', GetUnicodeChar('mu')),'FontSize', 18)
%     ylabel('d''','FontSize',18) 
%     box off
%     x = block_struct(i).DetectionRates.Amplitude;
%     xq = linspace(x(1), x(end));
%     yq = sigfun(block_struct(i).DetectionRates.dPrime, xq);
% %     plot(xq, yq, 'Color', [.4 .4 .4])
%     [~,b] = min(abs(yq-dprime_threshold));
%     plot([0 xq(b) xq(b)], [dprime_threshold, dprime_threshold -1], 'Color', [.4 .4 .4], 'LineStyle','--')
% 
% end
% https://www.mathworks.com/matlabcentral/answers/525387-how-to-use-a1-exp-x-b1-c1-2-formula-and-how-to-define-a-b-and-c-coefficents

% %% creating sliding window 
% %create structure
% %first and last 500 trials
% bin = struct(); ii = 1;
% monkey_name = 'Pinot';
% bin_type = {1 2};
% 
% for t = 1:length(bin_type)
% 
%     bin(ii).Monkey = monkey_name;
% 
%     bin(ii).bin_num = bin_type{t};
% 
%     bin_one = bigtable(1:500,:);
% 
%     bin_two = bigtable(end-500:end,:);
% 
%     bin(1).ResponseTable = bin_one;
%     bin(2).ResponseTable = bin_two;
% 
%     [detection_table_1, coeff_table_1] = AnalyzeDetectionTable(bin_one);
%     [detection_table_2, coeff_table_2] = AnalyzeDetectionTable(bin_two);
% 
%     bin(1).DetectionTable = detection_table_1; 
%     bin(2).DetectionTable = detection_table_2;
%     bin(1).CoeffTable = coeff_table_1; 
%     bin(2).CoeffTable = coeff_table_2;
% 
%     ii = ii + 1;
% end
% 
% [bigdetectiontable, bigcoefftable] = AnalyzeDetectionTable(bigtable);
% %% plotting and getting the threshold points
% 
% dprime_threshold = 1.35;
% % c(1) = rate of change, c(2) = x-offset, c(3) = multiplier, c(4) = offset
% sigfun = @(c,x) (c(3) .* (1./(1 + exp(-c(1).*(x-c(2)))))) + c(4); 
% binsize = (1:2);
% for i = 1:length(bin)
%     figure;hold on
%     x = bin(i).DetectionTable.StimAmp;
%     y = bin(i).DetectionTable.dPrime;
%     scatter(x,y,50, [.4 .4 .4], "filled")
%     plot(x,y, 'Color', [.4 .4 .4], 'LineStyle', ':')
%     xq = linspace(x(1), x(end));
%     yq = sigfun(bin(i).CoeffTable.dPrime, xq);
%     plot(xq, yq, 'Color', [.4 .4 .4])
%     [~,b] = min(abs(yq-dprime_threshold));
%     plot([0 xq(b) xq(b)], [dprime_threshold, dprime_threshold -1], 'Color', [.4 .4 .4], 'LineStyle','--')
%     text(x(end), 4-5*0.05, sprintf('d'' = %0.1f', xq(b)), 'VerticalAlignment', 'top', 'HorizontalAlignment', 'right')
%     set(gca, 'YLim', [-1 4], 'YTick', [-1:4])
%     ylabel('d'''); xlabel(sprintf('ICMS Amplitude (%sA)', GetUnicodeChar('mu')))
%     x0 = [100,-1];
%     ydata = [1 500];
%     something = lsqcurvefit(sigfun, x0, ydata, bin(2).CoeffTable);
% end
% %%
% figure; hold on
%  ax = gca;
%  ax.FontSize = 18;
% binsize = [500, 2000, 100];
% spacing = linspace(0,30, 2);
%  set(groot,'defaultLineMarkerSize',18);
% % luck = binsize(1);
% % unluck = bin(1).CoeffTable.dPrime(1,2);
% % scatter( bin(1).CoeffTable.dPrime(1,2),binsize(1))
% % scatter(bin(2).CoeffTable.dPrime(1,2),binsize