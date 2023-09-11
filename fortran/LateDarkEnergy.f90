module LateDE
    use DarkEnergyInterface
    use results
    use constants
    use classes
    implicit none

    private
    real(dl) :: grho_de_today

    type, extends(TDarkEnergyModel) :: TLateDE
        integer  :: model
        ! real(dl) :: winfty = -1.0_dl
        real(dl) :: w0, w1, w2, w3, w4, w5, w6, w7, w8, w9, w10
        real(dl) :: z1, z2, z3, z4, z5, z6, z7, z8, z9, z10 ! Binned w
        real(dl) :: fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8, fac9, fac10 ! Binned w factors

        contains        
        procedure :: w_de => TLateDE_w_de
        procedure :: grho_de => TLateDE_grho_de
        procedure :: Init => TLateDE_Init
        procedure :: ReadParams => TLateDE_ReadParams
        procedure :: PrintFeedback => TLateDE_PrintFeedback
        ! procedure :: BackgroundDensityAndPressure => TLateDE_density ! DHFS: Do I Need This ?
        procedure :: Effective_w_wa => TLateDE_Effective_w_wa   !VM: wont be called with CASARINI (our mod)
        procedure, nopass :: PythonClass => TLateDE_PythonClass
        procedure, nopass :: SelfPointer => TLateDE_SelfPointer
    end type TLateDE

    public TLateDE

    contains

    function TLateDE_w_de(this, a) result(w_de)
        class(TLateDE) :: this
        real(dl), intent(in) :: a    
        real(dl) :: w_de, z
        write(*,*)'1 model: ', this%model

        if (this%model == 1) then
            !'w_constant'
            w_de = this%w0
        else if (this%model == 2) then
            !'w0wa'
            w_de = this%w0 + this%w1*(1._dl - a)
        else if (this%model == 3) then
            !'3bins_w'
            z = 1._dl/a - 1._dl
            if (z < this%z1) then
                w_de = this%w0
            else if (z < this%z2) then
                w_de = this%w1
            else if (z < this%z3) then
                w_de = this%w2
            else
                w_de = this%w3
            end if    
        else if (this%model == 4) then
            !'5bins_w'
             z = 1._dl/a - 1._dl
            if (z < this%z1) then
                w_de = this%w0
            else if (z < this%z2) then
                w_de = this%w1
            else if (z < this%z3) then
                w_de = this%w2
            else if (z < this%z4) then
                w_de = this%w3
            else if (z < this%z5) then
                w_de = this%w4
            else
                w_de = this%w5
            end if     
        else if (this%model == 5) then
            !'10bins_w'
             z = 1._dl/a - 1._dl
            if (z < this%z1) then
                w_de = this%w0
            else if (z < this%z2) then
                w_de = this%w1
            else if (z < this%z3) then
                w_de = this%w2
            else if (z < this%z4) then
                w_de = this%w3
            else if (z < this%z5) then
                w_de = this%w4
            else if (z < this%z6) then
                w_de = this%w5
            else if (z < this%z7) then
                w_de = this%w6
            else if (z < this%z8) then
                w_de = this%w7
            else if (z < this%z9) then
                w_de = this%w8
            else if (z < this%z10) then
                w_de  = this%w9  
            else
                w_de = this%w10  
            end if    
        else
            stop "[Late Fluid DE] Invalid Dark Energy Model"   
        end if
    end function TLateDE_w_de

    function TLateDE_grho_de(this, a) result(grho_de)
        class(TLateDE) :: this
        real(dl), intent(in) :: a
        real(dl) :: grho_de, z, grho_de_today    
        write(*,*)'2 model: ', this%model

        ! Returns 8*pi*G * rho_de, no factor of a^4
        grho_de = 0
        if (this%model == 1) then
            !w = constant model
            grho_de = grho_de_today * a**(-3 * (1 + this%w0))
        else if (this%model == 2) then
            ! w0-wa model
            grho_de = grho_de_today * a**(3 * (1 + this%w0 + this%w1)) * exp(-3 * this%w1 * (1 - a))
        else if (this%model == 3) then
            ! Binned w model: 3 bins
            z = 1._dl/a - 1
            if (z < this%z1) then
                grho_de = grho_de_today * a**(-3 * (1 + this%w0))
            else if (a < this%z2) then
                grho_de = grho_de_today * this%fac1 * a**(-3 * (1 + this%w1))
            else if (z < this%z3) then
                grho_de = grho_de_today * this%fac2 * a**(-3 * (1 + this%w2))
            else
                grho_de = grho_de_today * this%fac3 * a**(-3 * (1 + this%z3))
            end if    
        else if (this%model == 4) then
            ! Binned w model: 5 bins
            z = 1._dl/a - 1
            if (z < this%z1) then
                grho_de = grho_de_today * a**(-3 * (1 + this%w0))
            else if (z < this%z2) then
                grho_de = grho_de_today * this%fac1 * a**(-3 * (1 + this%w1))
            else if (z < this%z3) then
                grho_de = grho_de_today * this%fac2 * a**(-3 * (1 + this%w2))
            else if (z < this%z4) then
                grho_de = grho_de_today * this%fac3 * a**(-3 * (1 + this%w3))
            else if (z < this%z5) then
                grho_de = grho_de_today * this%fac4 * a**(-3 * (1 + this%w4))
            else
                grho_de = grho_de_today * this%fac5 * a**(-3 * (1 + this%z5))
            end if    
        else if (this%model == 5) then
            ! Binned w model: 10 bins
            z = 1._dl/a - 1
            if (z < this%z1) then
                grho_de = grho_de_today * a**(-3 * (1 + this%w0))
            else if (z < this%z2) then
                grho_de = grho_de_today * this%fac1 * a**(-3 * (1 + this%w1))
            else if (z < this%z3) then
                grho_de = grho_de_today * this%fac2 * a**(-3 * (1 + this%w2))
            else if (z < this%z4) then
                grho_de = grho_de_today * this%fac3 * a**(-3 * (1 + this%w3))
            else if (z < this%z5) then
                grho_de = grho_de_today * this%fac4 * a**(-3 * (1 + this%w4))
            else if (z < this%z6) then
                grho_de = grho_de_today * this%fac5 * a**(-3 * (1 + this%w5))
            else if (z < this%z7) then
                grho_de = grho_de_today * this%fac6 * a**(-3 * (1 + this%w6))
            else if (z < this%z8) then
                grho_de = grho_de_today * this%fac7 * a**(-3 * (1 + this%w7))
            else if (z < this%z9) then
                grho_de = grho_de_today * this%fac8 * a**(-3 * (1 + this%w8)) 
            else if (z < this%z10) then
                grho_de = grho_de_today * this%fac9 * a**(-3 * (1 + this%w9))
            else
                grho_de = grho_de_today * this%fac10 * a**(-3 * (1 + this%z10))
            end if    
        else 
            stop "[Late Fluid DE] Invalid Dark Energy Model"
        end if
    end function TLateDE_grho_de

    subroutine TLateDE_Init(this, State)
        use classes
        use results
        class(TLateDE), intent(inout) :: this
        class(TCAMBdata), intent(in), target :: State
        write(*,*)'Model (It should be 1): ', this%model
        write(*,*)'----------------------------'
        write(*,*)'Eqn of State (w0, w1, ..., w10): '
        write(*,*) 'w0', this%w0
        write(*,*) 'w1', this%w1
        write(*,*) 'w2', this%w2
        write(*,*) 'w3', this%w3
        write(*,*) 'w4', this%w4
        write(*,*) 'w5', this%w5
        write(*,*) 'w6', this%w6
        write(*,*) 'w7', this%w7
        write(*,*) 'w8', this%w8
        write(*,*) 'w9', this%w9
        write(*,*) 'w10', this%w10
        write(*,*)'----------------------------'
        write(*,*)'Redshift (z1, z2 ,..., z10) : '
        write(*,*) 'z1', this%z1
        write(*,*) 'z2', this%z2
        write(*,*) 'z3', this%z3
        write(*,*) 'z4', this%z4
        write(*,*) 'z5', this%z5
        write(*,*) 'z6', this%z6
        write(*,*) 'z7', this%z7
        write(*,*) 'z8', this%z8
        write(*,*) 'z9', this%z9
        write(*,*) 'z10', this%z10

        ! this%is_cosmological_constant = .false.
        ! this%num_perturb_equations = 2 * 1
        
        ! if (this%model == 2 .or. this%model == 3 .or. this%model == 4 .or. this%model == 5) then
        !     this%num_perturb_equations = this%num_perturb_equations - 1 ! PPF only has one equation
        ! end if

        if (this%model == 3) then
            this%fac1 = (1+this%z1)**(3 * (this%w0 - this%w1))
            this%fac2 = this%fac1 * (1+this%z2)**(3 * (this%w1 - this%w2))
            this%fac3 = this%fac2 * (1+this%z3)**(3 * (this%w2 - this%w3))
        else if (this%model == 4) then
            this%fac1 = (1+this%z1)**(3 * (this%w0 - this%w1))
            this%fac2 = this%fac1 * (1+this%z2)**(3 * (this%w1 - this%w2))
            this%fac3 = this%fac2 * (1+this%z3)**(3 * (this%w2 - this%w3))
            this%fac4 = this%fac3 * (1+this%z4)**(3 * (this%w3 - this%w4))
            this%fac5 = this%fac4 * (1+this%z5)**(3 * (this%w4 - this%w5))
        else if (this%model == 5) then
            this%fac1 = (1+this%z1)**(3 * (this%w0 - this%w1))
            this%fac2 = this%fac1 * (1+this%z2)**(3 * (this%w1 - this%w2))
            this%fac3 = this%fac2 * (1+this%z3)**(3 * (this%w2 - this%w3))
            this%fac4 = this%fac3 * (1+this%z4)**(3 * (this%w3 - this%w4))
            this%fac5 = this%fac4 * (1+this%z5)**(3 * (this%w4 - this%w5))
            this%fac6 = this%fac5 * (1+this%z6)**(3 * (this%w5 - this%w6))
            this%fac7 = this%fac6 * (1+this%z7)**(3 * (this%w6 - this%w7))
            this%fac8 = this%fac7 * (1+this%z8)**(3 * (this%w7 - this%w8))
            this%fac9 = this%fac8 * (1+this%z9)**(3 * (this%w8 - this%w9))
            this%fac10 = this%fac9 * (1+this%z10)**(3 * (this%w9 - this%z10))
        end if  

        ! select type (State)
        !     type is (CAMBdata)
        !     grho_de_today = State%grhov
        ! end select        
    end subroutine TLateDE_Init

    subroutine TLateDE_ReadParams(this, Ini)
    !DHFS: copy from DarkEnergyInterface.f90 module - see TDarkEnergyEqnOfState_ReadParams
        use IniObjects
        use FileUtils
        class(TLateDE) :: this
        class(TIniFile), intent(in) :: Ini
        write(*,*)'4 model: ', this%model
    end subroutine TLateDE_ReadParams
    
    subroutine TLateDE_PrintFeedback(this, FeedbackLevel)
    !DHFS: copy from DarkEnergyInterface.f90 module - see TDarkEnergyEqnOfState_PrintFeedback
        class(TLateDE) :: this
        integer, intent(in) :: FeedbackLevel
        write(*,*)'5 model: ', this%model

        ! if (FeedbackLevel >0) write(*,'DHFS: what should I write here?')
    end subroutine TLateDE_PrintFeedback

    subroutine TLateDE_Effective_w_wa(this, w, wa)
        class(TLateDE), intent(inout) :: this
        real(dl), intent(out) :: w, wa
        write(*,*)'6 model: ', this%model

        if (this%model == 1) then
            !'w_constant'
            w  = this%w0
            wa = 0.0_dl
        else if (this%model == 2) then
            !'w0wa'
            w  = this%w0
            wa = this%w1
        else
            stop "[Late Fluid DE] Invalid Dark Energy Model (TLateDE_Effective_w_wa)"
        endif
    end subroutine TLateDE_Effective_w_wa

    function TLateDE_PythonClass()
        character(LEN=:), allocatable :: TLateDE_PythonClass

        TLateDE_PythonClass = 'LateDE'
    end function TLateDE_PythonClass

    subroutine TLateDE_SelfPointer(cptr,P)
        use iso_c_binding
        Type(c_ptr) :: cptr
        Type (TLateDE), pointer :: PType
        class (TPythonInterfacedClass), pointer :: P

        call c_f_pointer(cptr, PType)
        P => PType
    end subroutine TLateDE_SelfPointer

    ! subroutine TLateDE_density(this, grhov, a, grhov_t, w)
    !     ! Get grhov_t = 8*pi*G*rho_de*a**2 and (optionally) equation of state at scale factor a
    !     class(TLateDE), intent(inout) :: this
    !     real(dl), intent(in) :: grhov, a
    !     real(dl), intent(out) :: grhov_t
    !     real(dl), optional, intent(out) :: w
    !     write(*,*)'7 model: ', this%model

    !     if (a > 1e-10) then
    !         grhov_t = this%grho_de(a) * a**2
    !     else
    !         grhov_t = 0
    !     end if
    !     if (present(w)) then
    !         w = this%w_de(a)
    !     end if
    ! end subroutine TLateDE_density

end module LateDE