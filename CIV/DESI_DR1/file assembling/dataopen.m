%A = importdata('QSOdata.ssv',';');
%A = importdata('QSOdata.ssv',';',0);
B = importdata('QSOdata.csv',',');

% targetIDs = A.data(:,1);
% redshifts = A.data(:,2);
% wave = A.data(:,3);
% flux = A.data(:,4);
% ivar = A.data(:,5);
%% 
data = B{1};
rows = strsplit(data,'$');
n = numel(rows);
targetIDs = ones(1,n);
redshifts = ones(1,n);
for ii=1:n
    if isempty(rows{ii})
        continue
    end
    fields = strsplit(rows{ii}, ',');
    targetIDs(ii) = str2double(fields{1});
    redshifts(ii) = str2double(fields{2});
    wave{ii} = str2double(strsplit(fields{3}, ';'));    
    flux{ii} = str2double(strsplit(fields{4}, ';'));
    ivar{ii} = str2double(strsplit(fields{5}, ';'));
end
%% 
subplot(2,2,1)
plot(wave{5},flux{5})
subplot(2,2,2)
plot(wave{6},flux{6})
subplot(2,2,3)
plot(wave{7},flux{7})
subplot(2,2,4)
plot(wave{8},flux{8})
