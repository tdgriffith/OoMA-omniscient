function basis = gen_iH_basis(dim)
%% Overview
% Function to generate a basis for all self adjoint matricies of a given
% dim. Doesn't generate "Pauli like" structure, but should give identical
% results. Includes support for complex values

%diagonals
for i = 1:dim
    %%
    blank = zeros(dim,dim);
    blank(i,i) = 1j;
    diag{i} = blank;
end

%Off diagonal real
idx_start=2;
counter=0;
for i = 1:dim-1
    for ii = idx_start:dim
        if i ~= ii
            counter=counter+1;
            blank = zeros(dim,dim);
            blank(i,ii) = 1j;
            blank(ii,i) = 1j;
            off_diag_real{counter} = blank;
        end
    end
    idx_start=idx_start+1;
end

%Off diagonal complex
idx_start=2;
counter=0;
for i = 1:dim-1
    for ii = idx_start:dim
        if i ~= ii
            counter=counter+1;
            blank = zeros(dim,dim);
            blank(i,ii) = 1;
            blank(ii,i) = -1; %PLUS?
            off_diag_imag{counter} = blank;
        end
    end
    idx_start=idx_start+1;
end


basis=[diag,off_diag_real,off_diag_imag];



            
    