r_err = nan(3,8,10);
t_err = nan(3,8,10);
N=11;
n = 3;
for i = 1:N
    sequence = sprintf('%02d',i-1);
    [r_err(:,:,i),t_err(:,:,i)] = eval.evaluate_odometry(sequence,{'results_150_4_500','results_depth150_in4_ransac500_mono1'},{'ss','HX'});
end

r_avg = nan(n,N);
t_avg = nan(n,N);
for i=1:11
    for j = 1:n
        r = r_err(:,:,i);
        r = mean(r(j,~isnan(r(j,:))));
        
        t = t_err(:,:,i);
        t = mean(t(j,~isnan(t(j,:))));
        
        r_avg(j,i) = r;
        t_avg(j,i) = t;
    end
end

rr = mean(r_avg,2);
r_avg(:,12) = rr;

tt = mean(t_avg,2);
t_avg(:,12) = tt;

% This script runs two simple examples of latexTable.m
clc; clear input;

%% Example 1: using an numerical array as data input
fprintf('Example 1: using an array as data input\n\n');

% numeric values you want to tabulate:
% this field has to be an array or a MATLAB table
% in this example we use an array
input.data = r_avg;

% Optional fields:

% Set column labels (use empty string for no label):
input.tableColLabels = {'00','01','02','03','04','05','06','07','08','09','10','mean'};
% Set row labels (use empty string for no label):
input.tableRowLabels = {'SS','IO-S','IO-M'};

% Switch transposing/pivoting your table:
input.transposeTable = 0;

% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used

% Formatting-string to set the precision of the table values:
% For using different formats in different rows use a cell array like
% {myFormatString1,numberOfValues1,myFormatString2,numberOfValues2, ... }
% where myFormatString_ are formatting-strings and numberOfValues_ are the
% number of table columns or rows that the preceding formatting-string applies.
% Please make sure the sum of numberOfValues_ matches the number of columns or
% rows in input.tableData!
%
input.dataFormat = {'%.2e',12}; % three digits precision for first two columns, one digit for the last

% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = '-';

% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';

% Switch table borders on/off:
input.tableBorders = 1;

% LaTex table caption:
input.tableCaption = 'Rotation Error';

% LaTex table label:
input.tableLabel = 'fig:r_err';

% Switch to generate a complete LaTex document or just a table:
%input.makeCompleteLatexDocument = 1;

% call latexTable:
latex1 = latexTable(input);
fd = fopen('doc/eccv/r_err.tex','w+');
for i = 1:size(latex1,1)
    fprintf(fd,'%s\n',latex1{i,1});
end

% This script runs two simple examples of latexTable.m
clc; clear input;

%% Example 1: using an numerical array as data input
fprintf('Example 1: using an array as data input\n\n');

% numeric values you want to tabulate:
% this field has to be an array or a MATLAB table
% in this example we use an array
input.data = t_avg;

% Optional fields:

% Set column labels (use empty string for no label):
input.tableColLabels = {'00','01','02','03','04','05','06','07','08','09','10','mean'};
% Set row labels (use empty string for no label):
input.tableRowLabels = {'SS','IO-S','IO-M'};

% Switch transposing/pivoting your table:
input.transposeTable = 0;

% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used

% Formatting-string to set the precision of the table values:
% For using different formats in different rows use a cell array like
% {myFormatString1,numberOfValues1,myFormatString2,numberOfValues2, ... }
% where myFormatString_ are formatting-strings and numberOfValues_ are the
% number of table columns or rows that the preceding formatting-string applies.
% Please make sure the sum of numberOfValues_ matches the number of columns or
% rows in input.tableData!
%
input.dataFormat = {'%.3g',12}; % three digits precision for first two columns, one digit for the last

% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = '-';

% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';

% Switch table borders on/off:
input.tableBorders = 1;

% LaTex table caption:
input.tableCaption = 'Translation Error';

% LaTex table label:
input.tableLabel = 'fig:r_err';

% Switch to generate a complete LaTex document or just a table:
%input.makeCompleteLatexDocument = 1;

% call latexTable:
latex2 = latexTable(input);
fd = fopen('doc/eccv/t_err.tex','w+');
for i = 1:size(latex2,1)
    fprintf(fd,'%s\n',latex2{i,1});
end
