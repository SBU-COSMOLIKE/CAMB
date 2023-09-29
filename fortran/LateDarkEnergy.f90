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
        real(dl) :: w0,w1,w2,w3,w4,w5,w6,w7,w8,w9,w10
        real(dl) :: z1,z2,z3,z4,z5,z6,z7,z8,z9,z10
        real(dl) :: fac1,fac2,fac3,fac4,fac5,fac6,fac7,fac8,fac9,fac10
        contains
        procedure :: ReadParams => TLateDE_ReadParams
        procedure :: Init => TLateDE_Init
        procedure :: PrintFeedback => TLateDE_PrintFeedback
        procedure :: w_de => TLateDE_w_de
        procedure :: grho_de => TLateDE_grho_de
        procedure :: Effective_w_wa => TLateDE_Effective_w_wa   !VM: wont be called with CASARINI (our mod)         
        procedure, nopass :: SelfPointer => TLateDE_SelfPointer
        procedure :: BackgroundDensityAndPressure => TLateDE_density ! DHFS: Do I Need This ? If yes why, if not why
    end type TLateDE

    public TLateDE

    contains

    function TLateDE_w_de(this, a) result(w_de)
        class(TLateDE) :: this
        real(dl), intent(in) :: a    
        real(dl) :: w_de, z
        real(dl) :: w0,w1,w2,z0,z1,z2
        real(dl) :: wa0,waa0,wa1,waa1
        real(dl) :: Delta_z1, Delta_w1, Delta_z2, Delta_w2
        
        w_de = 0

        if (this%model == 1) then
            ! Constant w
            w_de = this%w0
        else if (this%model == 2) then
            ! CPL parametrization w0wa
            w_de = this%w0 + this%w1*(1._dl - a)
        else if (this%model == 3) then
            ! 3 bins w
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
            ! 5 bins w
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
            ! 10 bins w
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
        else if (this%model == 6) then
            ! 2 binned linear w(z). Linearity is in redshift, not in scale factor
            z = 1._dl/a - 1._dl
            if (z < this%z1) then
                w_de = this%w0 + (this%w1 - this%w0)/this%z1 * z
            else if (z < this%z2) then
                w_de = this%w1 + (this%w2 - this%w1)/(this%z2 - this%z1) * (z - this%z1)
            else
                w_de = this%w2
            end if  
        else if (this%model == 7) then
            ! 2 binned quadratic w(z)
            ! Definitions
            z0=0
            z1=this%z1
            z2=this%z2 
            w0=this%w0 
            w1=this%w1 
            w2=this%w2
            Delta_z1 = z1-z0
            Delta_w1 = w1-w0
            Delta_z2 = z2-z1
            Delta_w2 = w2-w1
            ! Boundary conditions -see my notes-
            wa0  = 2._dl*(Delta_w1/Delta_z1 - Delta_w2/Delta_z2)
            waa0 = -Delta_w1/Delta_z1**2 + 2._dl*Delta_w2/(Delta_z1*Delta_z2)
            wa1  =  2._dl*Delta_w2/Delta_z2
            waa1 =  -Delta_w2/Delta_z2**2
            ! Equation of state
            z = 1._dl/a - 1._dl
            if (z < this%z1) then
                w_de = w0 + wa0 * (z-z0) + waa0*(z-z0)**2
            else if (z < this%z2) then
                w_de = w1 + wa1 * (z-z1) + waa1*(z-z1)**2
            else
                w_de = w2
            end if                
        else
            stop "[Late Fluid DE] Invalid Dark Energy Model"   
        end if
    end function TLateDE_w_de

    function TLateDE_grho_de(this, a) result(grho_de)
        class(TLateDE) :: this
        real(dl), intent(in) :: a
        real(dl) :: grho_de, z    
        ! real(dl) :: wa0,wa1,wa2
        real(dl) :: alpha0,alpha1,alpha2

        real(dl) :: w0,w1,w2,w3,z0,z1,z2,z3
        real(dl) :: wa0,waa0,wa1,waa1,wa2,waa2
        real(dl) :: Delta_z1, Delta_w1, Delta_z2, Delta_w2, Delta_z3, Delta_w3  
        real(dl) :: A00,A10,A20, A01,A11,A21, A02,A12,A22      

        ! Definitions
        z0=0
        z1=this%z1
        z2=this%z2
        z3=this%z3 
        w0=this%w0 
        w1=this%w1 
        w2=this%w2
        w3=this%w3

        Delta_z1 = z1-z0
        Delta_z2 = z2-z1
        Delta_z3 = z3-z2
        Delta_w1 = w1-w0
        Delta_w2 = w2-w1
        Delta_w3 = w3-w2 

        ! Returns 8*pi*G * rho_de, no factor of a^4
        grho_de = 0

        if (this%model == 1) then
            ! w constant
            grho_de = grho_de_today * a**(-3 * (1 + this%w0))
        else if (this%model == 2) then
            ! CPL w0-wa
            grho_de = grho_de_today * a**(-3 * (1 + this%w0 + this%w1)) * exp(-3 * this%w1 * (1 - a))
        else if (this%model == 3) then
            ! 3 bins w
            z = 1._dl/a - 1
            if (z < this%z1) then
                grho_de = grho_de_today * a**(-3 * (1 + this%w0))
            else if (z < this%z2) then
                grho_de = grho_de_today * this%fac1 * a**(-3 * (1 + this%w1))
            else if (z < this%z3) then
                grho_de = grho_de_today * this%fac2 * a**(-3 * (1 + this%w2))
            else
                grho_de = grho_de_today * this%fac3 * a**(-3 * (1 + this%w3))
            end if    
        else if (this%model == 4) then
            ! 5 Bins w
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
                grho_de = grho_de_today * this%fac5 * a**(-3 * (1 + this%w5))
            end if    
        else if (this%model == 5) then
            ! 10 Bins w
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
                grho_de = grho_de_today * this%fac10 * a**(-3 * (1 + this%w10))
            end if    
        else if (this%model == 6) then
            ! 2 Bins linear w(z)           

            wa0 = Delta_w1/Delta_z1
            wa1 = Delta_w2/Delta_z2
            wa2 = Delta_w3/Delta_z3

            alpha0 = 3*(1+w0-wa0*(1+z0))      
            alpha1 = 3*(1+w1-wa1*(1+z1))
            alpha2 = 3*(1+w2-wa2*(1+z2))
            
            z = 1._dl/a - 1
            if (z < this%z1) then
                grho_de = grho_de_today * & 
                                        ((1+z)/(1+z0))**alpha0*exp(3*wa0*(z-z0))
            else if (z < this%z2) then
                grho_de = grho_de_today * &
                                        ((1+z1)/(1+z0))**alpha0*exp(3*wa0*(z1-z0)) * &
                                        ((1+z )/(1+z1))**alpha1*exp(3*wa1*(z -z1))
            else
                grho_de = grho_de_today * & 
                                        ((1+z1)/(1+z0))**alpha0*exp(3*wa0*(z1-z0)) * &
                                        ((1+z2)/(1+z1))**alpha1*exp(3*wa1*(z2-z1)) * &
                                        ((1+ z)/(1+z2))**alpha2*exp(3*wa2*(z -z2))
            end if 
        else if(this%model == 7) then 
            ! 2 bins quadratic w(z)
            z0=0
            z1=this%z1
            z2=this%z2 
            w0=this%w0 
            w1=this%w1 
            w2=this%w2
            Delta_z1 = z1-z0
            Delta_w1 = w1-w0
            Delta_z2 = z2-z1
            Delta_w2 = w2-w1   

            wa0  = 2._dl*(Delta_w1/Delta_z1 - Delta_w2/Delta_z2)
            waa0 = -Delta_w1/Delta_z1**2 + 2._dl*Delta_w2/(Delta_z1*Delta_z2)
            wa1  =  2._dl*Delta_w2/Delta_z2
            waa1 =  -Delta_w2/Delta_z2**2 
            wa2  = 0
            waa2 = 0

            A00 = 3._dl*(1+w0-wa0*z0+waa0*z0**2._dl)
            A10 = 3._dl*(wa0-2*waa0*z0)
            A20 = 3._dl*waa0

            A01 = 3._dl*(1+w1-wa1*z1+waa1*z1**2._dl)
            A11 = 3._dl*(wa1-2*waa1*z1)
            A21 = 3*waa1

            A02 = 3._dl*(1+w2-wa2*z2+waa2*z2**2._dl)
            A12 = 3._dl*(wa2-2*waa2*z2)
            A22 = 3._dl*waa2

            z = 1._dl/a - 1
            if (z < this%z1) then
                grho_de = grho_de_today * &
                          (((1+ z)/(1+z0))**(A00-A10+A20))*exp((A10-A20)*(z -z0)+A20*(z**2-z0**2)/2)
            else if (z < this%z2) then
                grho_de = grho_de_today * &
                          (((1+z1)/(1+z0))**(A00-A10+A20))*exp((A10-A20)*(z1-z0)+A20*(z1**2-z0**2)/2) * &
                          (((1+z )/(1+z1))**(A01-A11+A21))*exp((A11-A21)*(z -z1)+A21*( z**2-z1**2)/2)
            else
                grho_de = grho_de_today * & 
                          (((1+z1)/(1+z0))**(A00-A10+A20))*exp((A10-A20)*(z1-z0)+A20*(z1**2-z0**2)/2) * &
                          (((1+z2)/(1+z1))**(A01-A11+A21))*exp((A11-A21)*(z2-z1)+A21*(z2**2-z1**2)/2) * &
                          (((1+z )/(1+z2))**(A02-A12+A22))*exp((A12-A22)*(z -z2)+A22*( z**2-z2**2)/2)          
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
            this%fac10 = this%fac9 * (1+this%z10)**(3 * (this%w9 - this%w10))
        ! if (this%model == 6) then
            ! this%fac1 = 
            ! this%fac2 = 
            ! this%fac3 = 
        end if  

        select type (State)
            type is (CAMBdata)
            grho_de_today = State%grhov
        end select        
    end subroutine TLateDE_Init

    subroutine TLateDE_ReadParams(this, Ini)
        use IniObjects
        use FileUtils
        class(TLateDE) :: this
        class(TIniFile), intent(in) :: Ini
    end subroutine TLateDE_ReadParams
    
    subroutine TLateDE_PrintFeedback(this, FeedbackLevel)
        class(TLateDE) :: this
        integer, intent(in) :: FeedbackLevel

        if (FeedbackLevel >0) write(*,'("(w0, wa) = (", f8.5,", ", f8.5, ")")') &
        &   this%w_lam, this%wa
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

    subroutine TLateDE_SelfPointer(cptr,P)
        use iso_c_binding
        Type(c_ptr) :: cptr
        Type (TLateDE), pointer :: PType
        class (TPythonInterfacedClass), pointer :: P

        call c_f_pointer(cptr, PType)
        P => PType
    end subroutine TLateDE_SelfPointer

    subroutine TLateDE_density(this, grhov, a, grhov_t, w)
        ! Get grhov_t = 8*pi*G*rho_de*a**2 and (optionally) equation of state at scale factor a
        class(TLateDE), intent(inout) :: this
        real(dl), intent(in) :: grhov, a
        real(dl), intent(out) :: grhov_t
        real(dl), optional, intent(out) :: w

        if (a > 1e-10) then
            grhov_t = this%grho_de(a) * a**2
        else
            grhov_t = 0
        end if
        if (present(w)) then
            w = this%w_de(a)
        end if
    end subroutine TLateDE_density

end module LateDE