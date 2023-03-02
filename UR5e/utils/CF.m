function coef = CF(V,x,ind)
    N = length(V);
    coef = sym(zeros(N,1));
    if ind == 1
        var = x(1);
        for i = 1 : N            
            [dummy,target] = coeffs(V(i),var);
            idx = find(target==var,1);
            if isempty(idx)
                coef(i) = 0;
            else
                coef(i) = dummy(idx); % dummy(1) is constant term
            end
        end
    elseif ind == 2
        var1 = x(1); var2 = x(2);
        for i = 1 : N
            [dummy1,target1] = coeffs(V(i),var1);
            idx1 = find(target1==var1,1);
            if isempty(idx1)
                coef(i) = 0;
            else
                [dummy2,target2] = coeffs(dummy1(idx1),var2);
                idx2 = find(target2==var2,1);
                if isempty(idx2)
                    coef(i) = 0;
                else
                    coef(i) = dummy2(idx2);
                end
            end
        end
    elseif ind == 3
        var1 = x(1); var2 = x(2); var3 = x(3);
        for i = 1 : N
            [dummy1,target1] = coeffs(V(i),var1);
            idx1 = find(target1==var1,1);
            if isempty(idx1)
                coef(i) = 0;
            else
                [dummy2,target2] = coeffs(dummy1(idx1),var2);
                idx2 = find(target2==var2*var3,1);
                if isempty(idx2)
                    coef(i) = 0;
                else
                    coef(i) = dummy2(idx2);
                end
            end
        end
    elseif ind == 4
        var1 = x(1); var2 = x(2); var3 = x(3);
        for i = 1 : N
            [dummy1,target1] = coeffs(V(i),var1);
            idx1 = find(target1==var1,1);
            if isempty(idx1)
                coef(i) = 0;
            else
                [dummy2,target2] = coeffs(dummy1(idx1),var2);
                idx2 = find(target2==var2,1);
                if isempty(idx2)
                    coef(i) = 0;
                else
                    [dummy3,target3] = coeffs(dummy2(idx2),var3);
                    idx3 = find(target3==var3,1);
                    if isempty(idx3)
                        coef(i) = 0;
                    else
                        coef(i) = dummy3(idx3);
                    end
                end
            end
        end
    end
end