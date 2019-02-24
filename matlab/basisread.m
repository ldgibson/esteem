function basissetdef = basisread(basissetstring)
% Returns a cell array of structured arrays containing
%   information on shell type, exponents, and
%   contraction coefficients.

elements = {'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne'};
basis_func_num_dict = containers.Map;
basis_func_num_dict('STO-3G') = [1, 2];
basis_func_num_dict('6-31G') = [2, 3];
basis_func_num_dict('6-311G') = [3, 4];
basis_func_num_dict('cc-pVDZ') = [5, 9];

% throw error if not specifying one of the four basis sets
if ~ismember(basissetstring, basis_func_num_dict.keys)
    error('Invalid basis set. Please choose from the following options:\n  %s', ...
          'STO-3G, 6-31G, 6-311G, cc-pVDZ');
end

num_basis_func_array = basis_func_num_dict(basissetstring);

basissetdef = {};

filename = strcat('basissets/', basissetstring, '.basis');
basissetfile = fopen(filename, 'r');

id = 1;
while true
    line = fgetl(basissetfile);
    if ~ischar(line)  % end of file
        break
    end
    if startsWith(line, '!')  % skip comment line
        continue
    end
    if startsWith(line, '****')  % skip separators
        continue
    end
    line = strsplit(strtrim(line)); % remove trailing white space
    if ismember(line(1), elements)
        element = line(1);
        
        % if current element is H or He, only use valence basis functions
        if strcmp(element, 'H') || strcmp(element, 'He')
            num_basis_func = num_basis_func_array(1); % number of valence basis functions
        else
            num_basis_func = num_basis_func_array(2); % number of total basis functions
        end
        basissetdef{id} = struct('shelltype', cell(num_basis_func, 1), 'exponents', [], 'coeffs', []);
        for i = 1:length(basissetdef{id})
            info_line = fgetl(basissetfile);
            info_line = strsplit(strtrim(info_line));
            
            shell_type = info_line(1); % S, SP, P, D, etc.
            num_primitives = str2double(info_line(2));
            prefactor = str2double(info_line(3)); % unused
            
            exponents = zeros(1, num_primitives);
            if strcmp(shell_type, 'SP')
                contraction_coeffs = zeros(2, num_primitives);
            else
                contraction_coeffs = zeros(1, num_primitives);
            end
            
            % collect parameters for each primitive Gaussian in this shell
            for j = 1:num_primitives
                exp_and_coeff = str2num(fgetl(basissetfile));
                exponents(j) = exp_and_coeff(1);
                if strcmp(shell_type, 'SP')
                    contraction_coeffs(:, j) = exp_and_coeff(2:3);
                else
                    contraction_coeffs(j) = exp_and_coeff(2);
                end
            end
            % fill the struct
            basissetdef{id}(i).shelltype = shell_type;
            basissetdef{id}(i).exponents = exponents;
            basissetdef{id}(i).coeffs = contraction_coeffs;
        end
        id = id + 1;
    end
end

fclose(basissetfile);