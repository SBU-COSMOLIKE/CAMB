module LateDE
    use DarkEnergyInterface
    use results
    use constants
    use classes
    implicit none

    type, extends(TDarkEnergyModel) :: TLateDE
        real(dl) :: w_lam = -1_dl !VM not be used in Casarini except to init the search for effective constant w
        real(dl) :: wa    = 0._dl !VM not be used in Casarini except to init the search for effective constant w
        
        logical  :: no_perturbations = .false.
        integer  :: model
        real(dl) :: cs2_lam = 1_dl
        real(dl) :: winfty  = -1.0_dl
        real(dl) :: w0, w1, w2, w3, w4, w5, w6, w7, w8, w9
        real(dl) :: z1, z2, z3, z4, z5, z6, z7, z8, z9, z10 ! Binned w
        real(dl) :: fac1, fac2, fac3, fac4, fac5, fac6, fac7, fac8, fac9, fac10 ! Binned w factors

        contains        
        procedure :: w_de => TLateDE_w_de
        procedure :: grho_de => TLateDE_grho_de
        procedure :: Init => TLateDE_Init
        procedure :: ReadParams => TLateDE_ReadParams
        procedure :: PrintFeedback => TLateDE_PrintFeedback
        ! TODO - investigate if this function (Effective_w_wa) is used on Halofit Casarini 
        procedure :: Effective_w_wa => TLateDE_Effective_w_wa
    end type TLateDE

    contains

    function TLateDE_w_de(this, a) result(w_de)
        class(TLateDE) :: this
        real(dl), intent(in) :: a    
        real(dl) :: w_de, z

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
                w_de = this%winfty
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
                w_de = this%winfty
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
                w_de = this%winfty  
            end if    
        else
            stop "[Late Fluid DE] Invalid Dark Energy Model"   
        end if                              

    end function TLateDE_w_de

    function TLateDE_grho_de(this, a) result(grho_de)
    class(TLateDE) :: this
    real(dl), intent(in) :: a
    real(dl) :: grho_de, z, grho_de_today    

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
            grho_de = grho_de_today * this%fac3 * a**(-3 * (1 + this%winfty))
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
            grho_de = grho_de_today * this%fac5 * a**(-3 * (1 + this%winfty))
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
            grho_de = grho_de_today * this%fac10 * a**(-3 * (1 + this%winfty))
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
            this%fac10 = this%fac9 * (1+this%z10)**(3 * (this%w9 - this%winfty))
        end if    
    end subroutine TLateDE_Init

    subroutine TLateDE_ReadParams(this, Ini)
    !DHFS: copy from DarkEnergyInterface.f90 module - see TDarkEnergyEqnOfState_ReadParams
        use IniObjects
        use FileUtils
        class(TLateDE) :: this
        class(TIniFile), intent(in) :: Ini
    end subroutine TLateDE_ReadParams
    
    subroutine TLateDE_PrintFeedback(this, FeedbackLevel)
    !DHFS: copy from DarkEnergyInterface.f90 module - see TDarkEnergyEqnOfState_PrintFeedback
        class(TLateDE) :: this
        integer, intent(in) :: FeedbackLevel

        ! if (FeedbackLevel >0) write(*,'DHFS: what should I write here?')
    end subroutine TLateDE_PrintFeedback

    subroutine TLateDE_Effective_w_wa(this, w, wa)
        class(TLateDE), intent(inout) :: this
        real(dl), intent(out) :: w, wa

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

end module LateDE





!------------------------------------------------------------------



    ! function TLateDE_PythonClass()
    !     character(LEN=:), allocatable :: TLateDE_PythonClass

    !     TLateDE_PythonClass = 'LateDE'
    ! end function TLateDE_PythonClass

    ! subroutine TLateDE_SelfPointer(cptr,P)
    !     use iso_c_binding
    !     Type(c_ptr) :: cptr
    !     Type (TLateDE), pointer :: PType
    !     class (TPythonInterfacedClass), pointer :: P

    !     call c_f_pointer(cptr, PType)
    !     P => PType
    ! end subroutine TLateDE_SelfPointer

        ! subroutine TLateDE_density(this, grhov, a, grhov_t, w)
    !     ! Get grhov_t = 8*pi*G*rho_de*a**2 and (optionally) equation of state at scale factor a
    !     class(TLateDE), intent(inout) :: this
    !     real(dl), intent(in) :: grhov, a
    !     real(dl), intent(out) :: grhov_t
    !     real(dl), optional, intent(out) :: w

    !     if (a > 1e-10) then
    !         grhov_t = this%grho_de(a) * a**2
    !     else
    !         grhov_t = 0
    !     end if
    !     if (present(w)) then
    !         w = this%w_de(a)
    !     end if
    ! end subroutine TLateDE_density