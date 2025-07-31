detection_threshhold = 0.85;
detection_flag = p_c4 >= detection_threshhold;

[a,b] = size(detection_flag);
detected_TIDs = cell(a,1);
detected_zs = nan(a,b)%[];


detected_idx = 0;
for ii=1:a          %Loop over all QSOs that were run

    num_absorbers = nnz(detection_flag(ii,:));  %Number of detected absorbers in a QSO

    if num_absorbers >=1        %If there were any absorbers detected

        detected_idx = detected_idx + 1;    %Index in catalog output
        targetID = all_QSO_ID_dr1{ii};         %Find the QSOs target id

        detected_TIDs{detected_idx} = targetID; %Write the found target ID to an array
        %detected_TIDs = [detected_TIDs, fprintf('%s',targetID)];


        % Convert targetID to a string for comparison
        targetID_str = targetID; % targetID is already a char array

        % Initialize an array to hold the filtered indices
        filtered_idx = [];

        % Loop through each element in the cell array
        for i = 1:length(all_QSO_ID_dr1)
            % Compare the content of each cell with targetID
            if strcmp(all_QSO_ID_dr1{i}, targetID_str)
                filtered_idx = i; % Append the index if there's a match
            end
        end

        %filtered_idx = find(all_QSO_ID_dr1==targetID);
        sliced_zs = map_z_c4L2(filtered_idx,:);
           %Add other values like REW, col density here

        sliced_absorbers = sliced_zs(detection_flag(ii,:));
        
        for jj=1:numel(sliced_absorbers)
            detected_zs(detected_idx, jj) = sliced_absorbers(jj); % Store the detected absorbers
        end
    end
end
filled_cells = cellfun(@(x) ~isempty(x), detected_TIDs);
num_filled_cells = sum(filled_cells(:));

detected_TIDs(num_filled_cells+1:end) = [];
detected_zs(num_filled_cells+1:end,:) = [];

detected_TIDs
detected_zs

A = cell(numel(detected_TIDs), 8);
A(:,1) = detected_TIDs;
A(:,2:8) = num2cell(detected_zs);

filename = "CIV_catalog.csv";
writecell(A,filename)