load('Exams_T.mat');
n = 100;
nexams = 100;
nsps = length(sps);
nrs = length(rs);
nas = length(as);
timelimit = 5;
for i_sp = 1:nsps
    for i_r = 1:nrs
        for i_a = 1:nas
            sp = sps(i_sp);
            r = rs(i_r);
            a = as(i_a);
            for i_exam = 1
                B = Exams{i_sp,i_r,i_a,i_exam}{1};
                W = Exams{i_sp,i_r,i_a,i_exam}{2};
                A = W*W';
                n = size(A,1);
                H_0 = rand(size(W));
                options_ANLS.timelimit = timelimit;
                options_ANLS.Hinit = H_0;
                
                max_E = 3; %外部循环次数
             [H_anls,err_anls, ts_anls] = ANLS(A, size(H_0,2), options_ANLS, max_E);
            [H_bcd,err_bcd,ts_bcd] = BCD(A,50000,H_0,timelimit,max_E);
            end
        end
    end
end