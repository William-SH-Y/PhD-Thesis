function [n0, cyc, is_sj_crk, sj_strength, is_shear] = openfile_crk2(fname)
    filename = fopen(fname);
    C = textscan(filename,'%d %f %f %f %s (%f , %f , %f ) %f');
    cyc = cell2mat(C(1, 1));
    x = cell2mat(C(1, 2));
    y = cell2mat(C(1, 3));
    z = cell2mat(C(1, 4));
    crack_type = C{1, 5};
    nx = cell2mat(C(1, 6));
    ny = cell2mat(C(1, 7));
    nz = cell2mat(C(1, 8));
    sj_s = cell2mat(C(1, 9));

    temp = zeros(length(x), 1); %temp stores all the indices of SJ cracks from the origininal text file
                                %the FJ cracks will be 0
    is_sj_crk = zeros(length(x), 1);
    
    j = 0;
    for i = 1:length(x)
        if strcmpi(crack_type{i}, 'sj_Shear') || strcmpi(crack_type{i}, 'sj_Tension')
            j = j+1;
            temp(j) = i;
            is_sj_crk(i) = i;
        end
    end
    n0_n = nnz(temp); %nnz number of nonzero elements; n0_n nonzero number
    n0_ind = temp(1:n0_n); %n0_ind nonzero elements indices
    n0 = zeros(n0_n, 7); %n0 nonzero, as in SJ cracks only
    sj_strength = zeros(n0_n, 1);
    is_shear = zeros(n0_n, 1);
    for i = 1:n0_n
        n0(i, 1) = cyc(n0_ind(i)); %cycle which the SJ crack occurs
        n0(i, 2) = x(n0_ind(i));
        n0(i, 3) = y(n0_ind(i));
        n0(i, 4) = z(n0_ind(i));
        n0(i, 5) = nx(n0_ind(i));
        n0(i, 6) = ny(n0_ind(i));
        n0(i, 7) = nz(n0_ind(i));
        sj_strength(i) = sj_s(n0_ind(i));
        
        if strcmpi(crack_type{n0_ind(i)}, 'sj_Shear')
            is_shear(i) = 1;
        end
    end  
end