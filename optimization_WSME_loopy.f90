program optimization_WSME_loopy
    use defreal
    include 'nlopt.f'
    external dist
    real (kind=db)::parv(7),y(7)
    double precision d, params3(3), grad3(3),params5(5),grad5(5),minf_opt,eps,ds,dC,a,b
    double precision eps_i,ds_i,dC_i
    integer :: npar, flaggrad, ires, maxeval
    integer*8 opt,l
    REAL time_begin, time_end, tol1
    integer :: n_exp, flagpar, eval !n_exp is the number of experimental points, flagpar a flag to indicate if optimice 3 parameters or 5 (+a and b)
    real(kind=db), allocatable :: T_exp(:),C_exp(:) !an array to carry the experimental datapoints
    character(len=80):: expfile, simfile !the name of the .txt with the experimental datapoints and name of the .txt with the simulated data

    read(*,*) n_exp
    read(*,*) flagpar
    read(*,*) expfile
    read(*,*) simfile
    read(*,*) (parv(i),i=1,5)

    allocate(T_exp(1:n_exp),C_exp(1:n_exp))
    open(2,"Input/"//trim(expfile))
    do l=1,n_exp
        read(2,*) T_exp(l),C_exp(l)
    enddo
    close(2)
    open(20,"Output/optimizatoin_parameters.txt")

    eps=parv(1)
    ds=parv(2)
    dC=parv(3)
    a=parv(6)
    b=parv(7)

    flaggrad = 0 !flaggrad=0 (para no usar gradiente)
    grad = 0._db
    maxeval = 1000 !nº máximo de iteraciones permitidas
    tol1=1E-1
    opt = 0 
    eval=0
    CALL CPU_TIME ( time_begin ) 

    selectcase(flagpar)
    case(3)
        npar = 3 !nº parámetros a optimizar
        params3(1) = parv(1) !epsilon
        params3(2) = parv(2) !s
        params3(3) = parv(5) !deltaCp

        !parámetros que no se optimizan
        y(1)=parv(3) !epsilon_eff
        y(2)=parv(4) !I_solv
        y(3)=parv(6) !a
        y(4)=parv(7) !b
        
        write(*,*) 'Optimizando eps,s,dC' 

        call dist(d,n_exp,t_exp,C_exp,simfile,npar,params3,y,eval) !calcula distancia entre curva experimental y teorica
        write(*,*)  'd_min, dolddef iniciales = ',d,d*sqrt(real(n_exp)-1)/real(n_exp) 
            
        call nlo_create(opt, NLOPT_LN_NELDERMEAD, npar)
        call nlo_set_min_objective(ires, opt, dist, y)
        !call nlo_set_xtol_abs1(ires, opt, tol1)
        !call nlo_get_ftol_abs(tol1, opt)
        call nlo_set_maxeval(ires, opt, maxeval)
        call nlo_get_maxeval(maxeval, opt)
        call nlo_optimize(ires, opt, params3, minf_opt)

        if (ires.lt.0) then
            write(*,*) 'nlopt failed!'
        else
            write(*,*) 'Found min at eps,s,dC=', params3(1), params3(2), params3(3)
            write(*,*) 'd_min final, dolddef final = ',minf_opt, minf_opt*sqrt(real(n_exp)-1)/real(n_exp) 
            write(*,*) 'Number of iterations= ', eval-1
            !write(20,*) minf_opt,minf_opt*sqrt(real(n_exp)-1)/real(n_exp) , eval-1, params3(1), params3(2), params3(3),&
            !            &" !minf_opt, adjusted minf_opt, nº of evalutions, parameters to optimize"
            write(20,*) params3(1), params3(2),parv(3),parv(4), params3(3),parv(6),parv(7)," !all parameters in order"
            flush(20)
        endif
        call nlo_destroy(opt)

    case(5)
        npar = 5 !nº parámetros a optimizar
        params5(1) = parv(1) !epsilon
        params5(2) = parv(2) !s
        params5(3) = parv(5) !deltaCp
        params5(4) = parv(6) !a
        params5(5) = parv(7) !b

        !parámetros que no se optimizan
        y(1)=parv(3) !epsilon_eff
        y(2)=parv(4) !I_solv

        call dist(d,n_exp,T_exp,C_exp,simfile,npar,params5,y) !calcula distancia entre curva experimental y teorica
        write(*,*) 'd_min, dolddef iniciales = ',d,d*sqrt(real(n_exp)-1)/real(n_exp) 

        call nlo_create(opt, NLOPT_LN_NELDERMEAD, npar)
        call nlo_set_min_objective(ires, opt, dist, y)
        call nlo_set_maxeval(ires, opt, maxeval)
        call nlo_get_maxeval(maxeval, opt)
        call nlo_optimize(ires, opt, params5, minf_opt)

        if (ires.lt.0) then
           write(*,*) 'nlopt failed!'
        else
           write(*,*) 'found min at eps,s,dC=', params5(1), params5(2), params5(3)
           write(*,*) 'a,b=',params5(4),params5(5)
           write(*,*) 'd_min final = ', minf_opt
           write(*,*) 'nº iteraciones = ', eval-1
           !write(20,*) minf_opt,minf_opt*sqrt(real(n_exp)-1)/real(n_exp) , eval-1, params5(1), params5(2), params4(3), params5(4), params5(5)&
           !            &" !minf_opt, adjusted minf_opt, nº of evalutions, parameters to optimize"
           write(20,*) params5(1), params5(2), parv(3), parv(4), params5(3),params5(4),params5(5), " !all parameters in order"
           flush(20)
        endif
        call nlo_destroy(opt)

    end select

    CALL CPU_TIME ( time_end )
    WRITE (*,*) 'Time of operation was ', (time_end - time_begin)/60, ' minutes'
    write(*,*) '***** THE END *****'
    l=l+1 !don't know the function of l, maybe just a test
    

close(20)

end program optimization_WSME_loopy

subroutine dist(d,n_exp,T_exp, C_exp,simfile,npar,params,y,eval)
    use  defreal
    implicit none

    double precision :: d,params
    integer :: n_exp, npar, p, eval
    real (kind=db)::parv(7),y(7)
    real :: T_exp(n_exp), C_exp(n_exp)
    real :: T_sim(n_exp), C_sim(n_exp) !simulated data that WSME_genTden_loopy generates
    character(len=80) :: simfile
    character(len=200) :: command



    if(npar==3) then
        parv(1) = params(1) !eps
        parv(2) = params(2) !s
        parv(3)= y(1) !epsilon_eff
        parv(4) = y(2) !I_solv
        parv(5)= params(3) !deltaCp
        parv(6)= y(3) !a
        parv(7)= y(4) !b
    endif
  
    if(npar==5) then
        parv(1) = params(1) !eps
        parv(2) = params(2) !s
        parv(3)= y(1) !epsilon_eff
        parv(4) = y(2) !I_solv
        parv(5)= params(3) !deltaCp
        parv(6)= params(4) !a
        parv(7)= params(5) !b
    endif

    !Generate input file for WSME_genTden_loopy
    open(22,file="Input/opt_input_WSME.in")

    write(22,*) 128
    write(22,*) 1711
    write(22,*) 14400.0
    write(22,*) (parv(i), i=1,7)
    write(22,*) !temperaturas   HAY QUE ARREGLAR ESTO
    write(22,*) "0,0,1"
    write(22,*) 0 !number of (s,t) intervals
    write(22,*) "rCalpha.txt"
    write(22,*) "1DPX.map"
    write(22,*) "1DPX_elec.map"
    write(22,*) "cmapASA_1dpxmod_uf_HP.dat"
    write(22,*) "ct"
    write(22,*) ".false."
    write(22,*) ".true."
    write(22,*) ".false."
    write(22,*) ".false."
    write(22,*) ".false."
    write(22,*) ".false."
    write(22,*) ".false."
    write(22,*) ".false."
    write(22,*) ".false."
    write(22,*) ".false."
    write(22,*) ".true." !onlyCp
    write(22,*) ".false."
    write(22,*) ".false." !if .false. removes output from WSME_genTden in CMD


    call system('WSME_genTden_loopy.exe < Input/opt_input_WSME.in')
    open(94,file="Output/simCp.txt")
    read(94,*) (C_sim(p),p=1,n_exp)


    d=0.0_db
    do p=1,n_exp
        d=d+ABS(C_sim(p) - C_exp(p))*ABS(C_sim(p) - C_exp(p))
    end do 

    write(*,*) "Iteration number ",eval, " with distance ",d
    eval=eval+1
    
    close(22)
    close(94)
end subroutine