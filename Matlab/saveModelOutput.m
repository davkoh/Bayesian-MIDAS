function saveModelOutput(output, ...
                         Midas_type, gigg_type, sv_obs_type, ...
                         trend_type, midas_transformation_type)

    
    % Combine model name parts
    modelname = strcat('TREND_',trend_type, '_OBS_', sv_obs_type, '_PRIOR_', ...
                       Midas_type, '_GIGGTYPE_', gigg_type, '_TRANSFORM_', midas_transformation_type);

    output.modelname = modelname;
                 
    % Create folder name based on model definition
    foldername = fullfile('Output', modelname);
    
    % Ensure the folder exists; if not, create it
    if ~exist(foldername, 'dir')
        mkdir(foldername);
    end
    
    % Save the output structure
    save(fullfile(foldername, strcat('results', '.mat')), 'output');
    
    fprintf('Model output saved in folder: %s, with filename: results.mat\n', foldername);
end
