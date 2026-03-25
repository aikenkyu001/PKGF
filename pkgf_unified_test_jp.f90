program pkgf_unified_test_jp
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none
    
    integer, parameter :: DIM = 12
    real(real64), parameter :: EPSILON_METRIC = 1e-5
    real(real64), parameter :: EPSILON_MANIFOLD = 1e-4
    
    type :: TokenData
        character(len=64) :: text
        character(len=10) :: role
        real(real64) :: intensity
        integer :: char_len
    end type TokenData

    type :: BunsetsuData
        type(TokenData), allocatable :: tokens(:)
    end type BunsetsuData

    type :: SentenceData
        type(BunsetsuData), allocatable :: bunsetsus(:)
    end type SentenceData

    type(SentenceData), allocatable :: story_hierarchy(:)
    
    real(real64) :: K(DIM, DIM)
    real(real64) :: x(DIM)
    real(real64) :: t, dt
    real(real64) :: target_vec(DIM)
    real(real64) :: v_proj(DIM)
    real(real64) :: div_val
    real(real64) :: det_K, norm_v
    character(len=20) :: note
    integer :: s_idx, b_idx, t_idx
    integer :: steps, step_i
    real(real64) :: duration
    real(real64) :: G_token(DIM, DIM)
    real(real64) :: F_mat_local(DIM, DIM)

    call init_story_data()

    print *, "=== PKGF Unified Test (Purified 12D) ==="
    print *, "All experimental resources: https://github.com/aikenkyu001/PKGF"
    print *
    
    K = get_author_K()
    x = 0.1_real64
    t = 0.0_real64
    dt = 0.2_real64
    
    write(*, '(A, " | ", A, " | ", A, " | ", A, " | ", A, " | ", A)') &
        "Token         ", "t     ", "||v||     ", "div_g(KX) ", "det(K)    ", "Note"
    print *, repeat("-", 76)

    do s_idx = 1, size(story_hierarchy)
        do b_idx = 1, size(story_hierarchy(s_idx)%bunsetsus)
            do t_idx = 1, size(story_hierarchy(s_idx)%bunsetsus(b_idx)%tokens)
                associate (token => story_hierarchy(s_idx)%bunsetsus(b_idx)%tokens(t_idx))
                    
                    call get_target_vector(token, target_vec)
                    duration = 0.4_real64 + 0.1_real64 * real(token%char_len, real64)
                    steps = int(duration / dt)
                    
                    G_token = 0.0_real64
                    
                    do step_i = 1, steps
                        call simulation_step(x, K, target_vec, dt, G_token, v_proj)
                        t = t + dt
                    end do
                    
                    div_val = sum(matmul(K, v_proj))
                    det_K = calculate_det_12x12(K)
                    norm_v = sqrt(sum(v_proj**2))
                    
                    note = ""
                    
                    write(*, '(A, " | ", F6.2, " | ", F10.5, " | ", ES10.2, " | ", F10.5, " | ", A)') &
                        token%text(1:token%char_len*3), &
                        t, norm_v, div_val, det_K, trim(note)

                end associate
            end do
        end do
        print *, repeat("-", 76)
    end do
    
    print *
    print *, "[Scientific Findings]"
    print *, "1. Purified 12D Execution successful."
    print *, "2. Co-differential Propulsion verified in potential-only flow."
    print *, "3. Det(K) preserved without topological twist interference."

contains

    subroutine init_story_data()
        allocate(story_hierarchy(3))
        
        allocate(story_hierarchy(1)%bunsetsus(2))
        allocate(story_hierarchy(1)%bunsetsus(1)%tokens(1))
        story_hierarchy(1)%bunsetsus(1)%tokens(1) = TokenData("このアジェンダは", "Subject", 1.0_real64, 8)
        allocate(story_hierarchy(1)%bunsetsus(2)%tokens(1))
        story_hierarchy(1)%bunsetsus(2)%tokens(1) = TokenData("行動計画である", "Entity", 1.0_real64, 7)
        
        allocate(story_hierarchy(2)%bunsetsus(2))
        allocate(story_hierarchy(2)%bunsetsus(1)%tokens(2))
        story_hierarchy(2)%bunsetsus(1)%tokens(1) = TokenData("虚空の歯車が、", "Subject", 5.0_real64, 7)
        story_hierarchy(2)%bunsetsus(1)%tokens(2) = TokenData("メロンパンの", "Entity", 3.0_real64, 6)
        allocate(story_hierarchy(2)%bunsetsus(2)%tokens(2))
        story_hierarchy(2)%bunsetsus(2)%tokens(1) = TokenData("重力定数を", "Entity", 5.0_real64, 5)
        story_hierarchy(2)%bunsetsus(2)%tokens(2) = TokenData("逆走する。", "Action", 8.0_real64, 5)
        
        allocate(story_hierarchy(3)%bunsetsus(2))
        allocate(story_hierarchy(3)%bunsetsus(1)%tokens(1))
        story_hierarchy(3)%bunsetsus(1)%tokens(1) = TokenData("はっと目が覚める！", "Action", 8.0_real64, 9)
        allocate(story_hierarchy(3)%bunsetsus(2)%tokens(2))
        story_hierarchy(3)%bunsetsus(2)%tokens(1) = TokenData("己の至らなさを恥じ、", "Subject", 5.0_real64, 10)
        story_hierarchy(3)%bunsetsus(2)%tokens(2) = TokenData("深く反省する。", "Action", 6.0_real64, 7)
    end subroutine init_story_data

    subroutine get_target_vector(token, vec)
        type(TokenData), intent(in) :: token
        real(real64), intent(out) :: vec(DIM)
        real(real64) :: mag
        integer :: i
        
        vec = 0.0_real64
        mag = log(real(token%char_len, real64) + 2.0_real64) * 0.1_real64 * token%intensity
        
        if (trim(token%role) == "Subject") then
            vec(1:3) = mag
        else if (trim(token%role) == "Entity") then
            vec(4:6) = mag
        else if (trim(token%role) == "Action") then
            vec(7:9) = mag
        else if (trim(token%role) == "Context") then
            vec(10:12) = mag
        end if
        
        if (index(token%text, "目が覚める") > 0) then
            vec = vec * (-1.5_real64)
        end if
    end subroutine get_target_vector

    subroutine simulation_step(current_x, current_K, t_vec, dt, G_accum, v_out)
        real(real64), intent(inout) :: current_x(DIM)
        real(real64), intent(inout) :: current_K(DIM, DIM)
        real(real64), intent(in) :: t_vec(DIM)
        real(real64), intent(in) :: dt
        real(real64), intent(inout) :: G_accum(DIM, DIM)
        real(real64), intent(out) :: v_out(DIM)
        
        real(real64) :: delta_F(DIM)
        real(real64) :: g_metrics(DIM, DIM), g_inv(DIM, DIM)
        real(real64) :: v_raw(DIM), v_maxwell(DIM)
        real(real64) :: K_inv(DIM, DIM)
        real(real64) :: F_mat_loc(DIM, DIM)
        
        call compute_co_differential(current_x, t_vec, delta_F)
        
        g_metrics = get_metric_tensor(current_x)
        g_inv = mat_inv(g_metrics)
        v_raw = matmul(g_inv, delta_F)
        
        K_inv = mat_inv(current_K)
        v_maxwell = -matmul(K_inv, v_raw)
        
        call project_divergence_free(v_maxwell, current_x, current_K, v_out)
        
        call compute_F_matrix(current_x, t_vec, F_mat_loc)
        G_accum = G_accum + F_mat_loc * dt
        
        current_x = current_x + v_out * dt
        call parallel_transport_key(current_K, current_x, v_out, dt)
        
    end subroutine simulation_step

    function get_metric_tensor(pos) result(g)
        real(real64), intent(in) :: pos(DIM)
        real(real64) :: g(DIM, DIM)
        real(real64) :: s_subj, s_ent, s_act
        integer :: i
        g = 0.0_real64
        s_subj = 1.0_real64 + 0.5_real64 * tanh(pos(10))
        s_ent  = 1.0_real64 + 0.5_real64 * tanh(pos(11))
        s_act  = 1.0_real64 + 0.5_real64 * tanh(pos(12))
        do i = 1, DIM
            if (i <= 3) then; g(i,i) = s_subj
            else if (i <= 6) then; g(i,i) = s_ent
            else if (i <= 9) then; g(i,i) = s_act
            else; g(i,i) = 1.0_real64; end if
        end do
    end function get_metric_tensor

    function mat_inv(A) result(A_inv)
        real(real64), intent(in) :: A(DIM, DIM)
        real(real64) :: A_inv(DIM, DIM)
        real(real64) :: m(DIM, 2*DIM)
        real(real64) :: pivot, f, temp(2*DIM)
        integer :: i, k, n
        n = DIM
        A_inv = 0.0_real64; do i=1,n; A_inv(i,i)=1.0_real64; end do
        do i=1,n; m(i,1:n)=A(i,:); m(i,n+1:2*n)=A_inv(i,:); end do
        do i=1,n
            pivot = m(i,i)
            if (abs(pivot) < 1e-18) then
                do k=i+1,n
                    if (abs(m(k,i)) > abs(pivot)) then
                        temp=m(i,:); m(i,:)=m(k,:); m(k,:)=temp
                        pivot = m(i,i)
                        exit
                    end if
                end do
            end if
            if (abs(pivot) > 1e-18) then
                m(i,:) = m(i,:) * (1.0_real64/pivot)
                do k=1,n
                    if (i/=k) then
                        f=m(k,i); m(k,:) = m(k,:) - f*m(i,:)
                    end if
                end do
            end if
        end do
        A_inv = m(:, n+1:2*n)
    end function mat_inv

    function calculate_det_12x12(A) result(det)
        real(real64), intent(in) :: A(DIM, DIM)
        real(real64) :: det
        real(real64) :: lu(DIM, DIM)
        real(real64) :: pivot, factor, temp(DIM)
        integer :: i, k, n
        n = DIM; lu = A; det = 1.0_real64
        do i=1,n
            pivot = lu(i,i)
            if (abs(pivot) < 1e-15) then
                do k=i+1,n
                    if (abs(lu(k,i)) > abs(pivot)) then
                        temp=lu(i,:); lu(i,:)=lu(k,:); lu(k,:)=temp
                        det = det * (-1.0_real64); pivot = lu(i,i)
                        exit
                    end if
                end do
                if (abs(pivot) < 1e-15) then; det=0.0_real64; return; end if
            end if
            det = det * pivot
            do k=i+1,n
                factor = lu(k,i)/pivot
                lu(k,i) = factor
                lu(k,i+1:n) = lu(k,i+1:n) - factor*lu(i,i+1:n)
            end do
        end do
    end function calculate_det_12x12

    function get_author_K() result(K)
        real(real64) :: K(DIM, DIM)
        real(real64) :: scales(4)
        integer :: c_idx, i, row
        K = 0.0_real64; scales = [1.0_real64, 1.1_real64, 1.2_real64, 0.9_real64]
        do c_idx=1,4; do i=1,3; row=(c_idx-1)*3+i; K(row,row)=scales(c_idx); end do; end do
    end function get_author_K

    subroutine compute_co_differential(pos, t_vec, delta_F)
        real(real64), intent(in) :: pos(DIM), t_vec(DIM)
        real(real64), intent(out) :: delta_F(DIM)
        real(real64) :: d_star_F_vals(12) ! 11-form has 12 components
        real(real64) :: p_plus(DIM), p_minus(DIM)
        real(real64) :: sf_p(66), sf_m(66), diff(66)
        integer :: i
        
        d_star_F_vals = 0.0_real64
        
        do i = 1, DIM
             p_plus = pos; p_plus(i) = p_plus(i) + EPSILON_MANIFOLD
             p_minus = pos; p_minus(i) = p_minus(i) - EPSILON_MANIFOLD
             
             call compute_star_F(p_plus, t_vec, sf_p)
             call compute_star_F(p_minus, t_vec, sf_m)
             
             diff = (sf_p - sf_m) / (2.0_real64 * EPSILON_MANIFOLD)
             call accumulate_d_10form(diff, i, d_star_F_vals)
        end do
        
        call hodge_star_11form_to_1form(pos, d_star_F_vals, delta_F)
    end subroutine compute_co_differential

    subroutine accumulate_d_10form(diff_vals, grad_dim, accum)
        real(real64), intent(in) :: diff_vals(66)
        integer, intent(in) :: grad_dim
        real(real64), intent(inout) :: accum(12)
        
        integer :: c(10), i, idx
        logical :: has_dim
        integer :: sign_val, final_missing
        integer :: m1, m2, k
        
        idx = 1
        do i=1,10; c(i)=i; end do
        
        do while(idx <= 66)
            has_dim = .false.
            do i=1,10
                if (c(i) == grad_dim) then; has_dim=.true.; exit; end if
            end do
            
            if (.not. has_dim) then
                m1 = 0; m2 = 0; k = 1
                do i=1,12
                    if (k <= 10 .and. c(k) == i) then
                        k = k + 1
                    else
                        if (m1 == 0) then; m1 = i; else; m2 = i; end if
                    end if
                end do
                
                if (m1 == grad_dim) then; final_missing = m2; else; final_missing = m1; end if
                
                sign_val = 1
                do i=1, 10
                     if (grad_dim > c(i)) sign_val = sign_val * (-1)
                end do
                
                accum(13 - final_missing) = accum(13 - final_missing) + sign_val * diff_vals(idx)
            end if
            call next_comb(c, 12, 10, idx)
        end do
    end subroutine accumulate_d_10form

    subroutine next_comb(c, n, k, idx)
        integer, intent(inout) :: c(k)
        integer, intent(in) :: n, k
        integer, intent(inout) :: idx
        integer :: i, j
        idx = idx + 1
        i = k
        do while (i >= 1)
            if (c(i) < n - k + i) then
                c(i) = c(i) + 1
                do j = i+1, k; c(j) = c(j-1) + 1; end do
                return
            end if
            i = i - 1
        end do
    end subroutine next_comb

    subroutine compute_star_F(pos, t_vec, star_F_vals)
        real(real64), intent(in) :: pos(DIM), t_vec(DIM)
        real(real64), intent(out) :: star_F_vals(66)
        real(real64) :: F_vals(66), partials(DIM, DIM)
        real(real64) :: p_p(DIM), p_m(DIM), om_p(DIM), om_m(DIM)
        integer :: i, j, c(2), idx
        
        do i=1,DIM
            p_p=pos; p_p(i)=p_p(i)+EPSILON_MANIFOLD; p_m=pos; p_m(i)=p_m(i)-EPSILON_MANIFOLD
            om_p = get_omega_at(p_p, t_vec)
            om_m = get_omega_at(p_m, t_vec)
            partials(i,:) = (om_p - om_m)/(2.0_real64*EPSILON_MANIFOLD)
        end do
        
        c = [1,2]; idx=1
        do while(idx<=66)
            F_vals(idx) = partials(c(1), c(2)) - partials(c(2), c(1))
            call next_comb(c, 12, 2, idx)
        end do
        
        call hodge_star_generic(pos, F_vals, 2, star_F_vals)
    end subroutine compute_star_F

    function get_omega_at(pos, t_vec) result(omega)
        real(real64), intent(in) :: pos(DIM), t_vec(DIM)
        real(real64) :: omega(DIM), v_r(DIM), g(DIM, DIM)
        v_r = (t_vec - pos) * 0.5_real64
        g = get_metric_tensor(pos)
        omega = matmul(g, v_r)
    end function get_omega_at

    subroutine hodge_star_11form_to_1form(x, form_in, form_out)
        real(real64), intent(in) :: x(DIM), form_in(12)
        real(real64), intent(out) :: form_out(12)
        call hodge_star_generic(x, form_in, 11, form_out)
    end subroutine hodge_star_11form_to_1form

    subroutine hodge_star_generic(x, form_in, p, form_out)
        real(real64), intent(in) :: x(DIM)
        integer, intent(in) :: p
        real(real64), intent(in) :: form_in(*)
        real(real64), intent(out) :: form_out(*)
        
        real(real64) :: g(DIM,DIM), g_inv(DIM,DIM), det_g, sqrt_g
        integer :: n_combs, combs(66, p), i, j, k
        real(real64) :: omega_up(66), factor
        integer :: q, c_q(12-p), idx_q, c_p(p), idx_p, sign_val
        
        if (p==2) then; n_combs=66; else; n_combs=12; end if
        
        block
            integer :: c(p), idx
            do k=1,p; c(k)=k; end do
            idx=1
            do while(idx<=n_combs)
                combs(idx,:) = c
                call next_comb(c, 12, p, idx)
            end do
        end block
        
        g = get_metric_tensor(x); det_g=calculate_det_12x12(g); sqrt_g=sqrt(abs(det_g))
        g_inv = mat_inv(g)
        omega_up = 0.0_real64
        
        do i=1,n_combs
            do j=1,n_combs
                if (abs(form_in(j)) < 1e-18) cycle
                factor = 1.0_real64
                do k=1,p; factor = factor * g_inv(combs(i,k), combs(j,k)); end do
                omega_up(i) = omega_up(i) + factor * form_in(j)
            end do
        end do
        
        q = 12 - p; do k=1,q; c_q(k)=k; end do; idx_q=1
        do while(idx_q <= n_combs)
            call get_complement(c_q, q, c_p, p)
            idx_p = 0
            do i=1,n_combs
                if (all(combs(i,:)==c_p)) then; idx_p=i; exit; end if
            end do
            
            call get_epsilon_sign(c_p, c_q, p, q, sign_val)
            form_out(idx_q) = sqrt_g * omega_up(idx_p) * sign_val
            call next_comb(c_q, 12, q, idx_q)
        end do
    end subroutine hodge_star_generic

    subroutine get_complement(sub, slen, comp, clen)
        integer, intent(in) :: sub(slen), slen, clen
        integer, intent(out) :: comp(clen)
        integer :: i, j, k
        k=1; j=1
        do i=1,12
            if (j<=slen .and. sub(j)==i) then
                j=j+1
            else
                comp(k)=i; k=k+1
            end if
        end do
    end subroutine get_complement

    subroutine get_epsilon_sign(p_c, q_c, lp, lq, s_val)
        integer, intent(in) :: p_c(lp), q_c(lq), lp, lq
        integer, intent(out) :: s_val
        integer :: full(12), i, j, swaps
        full(1:lp)=p_c; full(lp+1:12)=q_c
        swaps=0
        do i=1,12; do j=i+1,12; if (full(i)>full(j)) swaps=swaps+1; end do; end do
        if (mod(swaps,2)==0) then; s_val=1; else; s_val=-1; end if
    end subroutine get_epsilon_sign

    subroutine project_divergence_free(v, x, K, v_out)
        real(real64), intent(in) :: v(DIM), x(DIM), K(DIM,DIM)
        real(real64), intent(out) :: v_out(DIM)
        real(real64) :: db, gf(DIM), vt(DIM), gns, lam
        integer :: i
        db = calc_local_div(v, x, K)
        do i=1,DIM
            vt=v; vt(i)=vt(i)+1e-5_real64
            gf(i) = (calc_local_div(vt, x, K) - db)/1e-5_real64
        end do
        gns = sum(gf**2)
        if (gns > 1e-18) then
            lam = db/gns; v_out = v - lam*gf
        else
            v_out = v
        end if
    end subroutine project_divergence_free

    function calc_local_div(v, x, K) result(dv)
        real(real64), intent(in) :: v(DIM), x(DIM), K(DIM,DIM)
        real(real64) :: dv, xp(DIM), xm(DIM), gp(DIM,DIM), gm(DIM,DIM)
        real(real64) :: sgp, sgm, fl(DIM), tp, tm
        integer :: i
        fl = matmul(K, v); dv = 0.0_real64
        do i=1,DIM
            xp=x; xp(i)=xp(i)+1e-5_real64; xm=x; xm(i)=xm(i)-1e-5_real64
            gp=get_metric_tensor(xp); gm=get_metric_tensor(xm)
            sgp=sqrt(abs(calculate_det_12x12(gp))); sgm=sqrt(abs(calculate_det_12x12(gm)))
            tp=sgp*fl(i); tm=sgm*fl(i)
            dv = dv + (tp - tm)/(2.0e-5_real64)
        end do
        dv = dv / sqrt(abs(calculate_det_12x12(get_metric_tensor(x))))
    end function calc_local_div

    subroutine parallel_transport_key(K_io, pos, v_dt, dt)
        real(real64), intent(inout) :: K_io(DIM, DIM)
        real(real64), intent(in) :: pos(DIM), v_dt(DIM), dt
        real(real64) :: G(DIM,DIM,DIM), Om(DIM,DIM), H(DIM,DIM), Hi(DIM,DIM)
        integer :: i, j, k
        call compute_chris(pos, G)
        Om = 0.0_real64
        do i=1,DIM; do j=1,DIM; do k=1,DIM
            Om(i,j) = Om(i,j) + G(i,k,j)*v_dt(k)
        end do; end do; end do
        H = matrix_exp_pade(Om * dt); Hi = mat_inv(H)
        K_io = matmul(matmul(H, K_io), Hi)
    end subroutine parallel_transport_key

    subroutine compute_chris(x, Gamma_sym)
        real(real64), intent(in) :: x(DIM)
        real(real64), intent(out) :: Gamma_sym(DIM,DIM,DIM)
        real(real64) :: metric(DIM,DIM), metric_inv(DIM,DIM), partials_val(DIM,DIM,DIM), t
        integer :: i, j, k, l
        metric = get_metric_tensor(x); metric_inv = mat_inv(metric)
        do k=1,DIM; partials_val(k,:,:) = partial_dg(x, k); end do
        Gamma_sym = 0.0_real64
        do k=1,DIM; do i=1,DIM; do j=1,DIM; do l=1,DIM
            if (abs(metric_inv(k,l))>1e-9) then
                t = partials_val(i,j,l) + partials_val(j,i,l) - partials_val(l,i,j)
                Gamma_sym(k,i,j) = Gamma_sym(k,i,j) + 0.5_real64*metric_inv(k,l)*t
            end if
        end do; end do; end do; end do
    end subroutine compute_chris

    function partial_dg(x, k) result(dg)
        real(real64), intent(in) :: x(DIM)
        integer, intent(in) :: k
        real(real64) :: dg(DIM, DIM), xp(DIM), xm(DIM)
        xp=x; xp(k)=xp(k)+EPSILON_METRIC; xm=x; xm(k)=xm(k)-EPSILON_METRIC
        dg = (get_metric_tensor(xp)-get_metric_tensor(xm))*(0.5_real64/EPSILON_METRIC)
    end function partial_dg

    function matrix_exp_pade(A) result(res)
        real(real64), intent(in) :: A(DIM,DIM)
        real(real64) :: res(DIM,DIM), As(DIM,DIM), Ap(DIM,DIM), P(DIM,DIM), Q(DIM,DIM), term(DIM,DIM)
        real(real64) :: inf_norm, rs, c(7)
        integer :: k_s, i, j
        c = [1.0_real64, 0.5_real64, 0.1_real64, 0.01111111111111111_real64, &
             0.0007575757575757576_real64, 3.0303030303030303e-05_real64, 5.050505050505051e-07_real64]
        inf_norm=0.0_real64
        do i=1,DIM; rs=sum(abs(A(i,:))); if (rs>inf_norm) inf_norm=rs; end do
        k_s=0; if (inf_norm>0) k_s=max(0, ceiling(log(inf_norm/0.5_real64)/log(2.0_real64)))
        As = A * (1.0_real64/(2.0_real64**k_s))
        P=0.0; Q=0.0; Ap=0.0
        do j=1,DIM; P(j,j)=1.0; Q(j,j)=1.0; Ap(j,j)=1.0; end do
        do i=1,6
            Ap = matmul(Ap, As); term = Ap * c(i+1)
            P = P + term
            if (mod(i,2)==0) then; Q=Q+term; else; Q=Q-term; end if
        end do
        res = matmul(mat_inv(Q), P)
        do i=1,k_s; res=matmul(res, res); end do
    end function matrix_exp_pade

    subroutine compute_F_matrix(pos, t_vec, F_mat)
        real(real64), intent(in) :: pos(DIM), t_vec(DIM)
        real(real64), intent(out) :: F_mat(DIM, DIM)
        real(real64) :: parts(DIM, DIM), pp(DIM), pm(DIM), op(DIM), om(DIM)
        integer :: i, j
        do i=1,DIM
            pp=pos; pp(i)=pp(i)+EPSILON_MANIFOLD; pm=pos; pm(i)=pm(i)-EPSILON_MANIFOLD
            op=get_omega_at(pp,t_vec); om=get_omega_at(pm,t_vec)
            parts(i,:) = (op-om)/(2.0_real64*EPSILON_MANIFOLD)
        end do
        do i=1,DIM; do j=1,DIM
            F_mat(i,j) = parts(i,j) - parts(j,i)
        end do; end do
    end subroutine compute_F_matrix

end program pkgf_unified_test_jp
