%% Least Squares a Given Matrix to the "Nearest Self-Adjoint Matrix"

function matrix_SA=least_squares_SA(A)
    [rows,cols]=size(A);
    
    if rows~=cols
        disp('Not square, check matrix dummy')
    end
    
    off_diag_terms=trianglenumbers(rows-1);
    diag_terms=rows;
    
    x_length=(diag_terms)+(2*off_diag_terms);
    B_length=(diag_terms)+(4*off_diag_terms);
    B=[];
    B_diag=[];
    B_offdiag=[];
    diag_locs=[];
    
    for i=1:rows
        for ii=1:rows
            if i==ii
                B_diag=[B_diag;real(A(i,ii))];
                diag_locs=[diag_locs,length(B_diag)]
            elseif i~=ii
                B_offdiag=[B_offdiag;real(A(i,ii))];
                B_offdiag=[B_offdiag;imag(A(i,ii))];
                B_offdiag=[B_offdiag;real(A(ii,i))];
                B_offdiag=[B_offdiag;imag(A(ii,i))];
            end
        end
    end
    

    for i=1:rows
        B_diag=[B_diag;real(A(i,i))];
        diag_locs=[diag_locs,length(B_diag)];
    end
    j=1;
    k=2;
    while j+k<2*rows
        if j~=k
            B_offdiag=[B_offdiag;real(A(j,k))];
            B_offdiag=[B_offdiag;imag(A(j,k))];
            B_offdiag=[B_offdiag;real(A(k,j))];
            B_offdiag=[B_offdiag;imag(A(k,j))];
        end
    end
    B=[B_diag;B_offdiag];
    clear i ii
    
    A_solve=zeros(B_length,x_length);
    block=[1 0;0 1;1 0;0 -1];
    A_solve(1:rows,1:rows)=eye(rows);
    j=rows+1;
    k=rows+1;
    for i=1:B_length/4
        A_solve(j:j+3,k:k+1)=block;
        j=j+4;
        k=k+2;
    end
    

            
            
            
            

        
    %working for 2x, not for 3x, something with block I think



end