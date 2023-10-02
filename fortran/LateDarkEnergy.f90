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
        real(dl) :: w0,w1,w2,w3,w4,w5,w6,w7,w8,w9
        real(dl) :: z1,z2,z3,z4,z5,z6,z7,z8,z9,z10
        real(dl) :: z0 = 0.0_dl
        contains
        procedure :: ReadParams => TLateDE_ReadParams
        procedure :: Init => TLateDE_Init
        procedure :: PrintFeedback => TLateDE_PrintFeedback
        procedure :: w_de => TLateDE_w_de
        procedure :: grho_de => TLateDE_grho_de
        procedure :: Effective_w_wa => TLateDE_Effective_w_wa         !VM: wont be called with CASARINI (our mod)         
        procedure, nopass :: SelfPointer => TLateDE_SelfPointer
        procedure :: BackgroundDensityAndPressure => TLateDE_density ! DHFS: Do I Need This ? If yes why, if not why
    end type TLateDE

    public TLateDE

    contains

    function TLateDE_w_de(this, a) result(w_de)
        class(TLateDE) :: this
        real(dl), intent(in) :: a    
        real(dl) :: w_de, z
        real(dl) :: wa0,waa0,wa1,waa1
        real(dl) :: Delta_z1, Delta_w1, Delta_z2, Delta_w2
        
        w_de = 0
        z = 1._dl/a - 1._dl

        if (this%model == 1) then
            ! Constant w
            w_de = this%w0
        else if (this%model == 2) then
            ! CPL parametrization w0wa
            w_de = this%w0 + this%w1*(1._dl - a)
        else if (this%model == 3) then
            ! 3 bins w
            if (z < this%z1) then
                w_de = this%w0
            else if (z < this%z2) then
                w_de = this%w1
            else if (z < this%z3) then
                w_de = this%w2
            else
                w_de = -1.0
            end if    
        else if (this%model == 4) then
            ! 5 bins w
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
                w_de = -1.0
            end if     
        else if (this%model == 5) then
            ! 10 bins w
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
                w_de = -1.0
            end if
        else if (this%model == 6) then
            !DHF 2 linear bins (in redshift) w(z).
            if (z < this%z1) then
                w_de = this%w0 + (this%w1 - this%w0)/this%z1 * z
            else if (z < this%z2) then
                w_de = this%w1 + ((-1.0) - this%w1)/(this%z2 - this%z1) * (z - this%z1)
            else
                w_de = (-1.0)
            end if  
        else if (this%model == 7) then
            ! 2 quadratic bins (in redshift)  w(z)
            Delta_z1 = this%z1 - 0.0
            Delta_w1 = this%w1 - this%w0
            Delta_z2 = this%z2 - this%z1
            Delta_w2 = (-1.0)  - this%w1
            
            !DHF: Boundary conditions -see my notes-
            wa0  = 2._dl*(Delta_w1/Delta_z1 - Delta_w2/Delta_z2)
            waa0 = -Delta_w1/Delta_z1**2 + 2._dl*Delta_w2/(Delta_z1*Delta_z2)
            wa1  =  2._dl*Delta_w2/Delta_z2
            waa1 =  -Delta_w2/Delta_z2**2
            
            if (z < this%z1) then
                w_de = this%w0 + wa0 * (z - 0.0) + waa0*(z - 0.0)**2
            else if (z < this%z2) then
                w_de = this%w1 + wa1 * (z - this%z1) + waa1*(z - this%z1)**2
            else
                w_de = -1.0
            end if                
        else
            stop "[Late Fluid DE @TLateDE_w_de] Invalid Dark Energy Model"   
        end if
    end function TLateDE_w_de

    function TLateDE_grho_de(this, a) result(grho_de)
        class(TLateDE) :: this
        real(dl), intent(in) :: a
        real(dl) :: grho_de, z    
        real(dl) :: alpha0, alpha1, alpha2
        real(dl) :: wa0,waa0,wa1,waa1,wa2,waa2
        real(dl) :: fac1,fac2,fac3,fac4,fac5,fac6,fac7,fac8,fac9,fac10
        real(dl) :: Delta_z1, Delta_w1, Delta_z2, Delta_w2, Delta_z3, Delta_w3  
        real(dl) :: A00,A10,A20, A01,A11,A21, A02,A12,A22      

        ! Returns 8*pi*G * rho_de, no factor of a^4
        grho_de = 0
        z = 1._dl/a - 1

        if (this%model == 1) then
            ! w constant
            grho_de = grho_de_today * a**(-3 * (1 + this%w0))
        else if (this%model == 2) then
            ! CPL w0-wa
            grho_de = grho_de_today * a**(-3 * (1 + this%w0 + this%w1)) * exp(-3 * this%w1 * (1 - a))
        else if (this%model == 3) then
            ! 3 bins w
            fac1 = (1.0+this%z1)**(3 * (this%w0 - this%w1))
            fac2 = fac1 * (1.0+this%z2)**(3 * (this%w1 - this%w2))
            fac3 = fac2 * (1.0+this%z3)**(3 * (this%w2 - (-1)))

            if (z < this%z1) then
                grho_de = grho_de_today * a**(-3 * (1 + this%w0))
            else if (z < this%z2) then
                grho_de = grho_de_today * fac1 * a**(-3 * (1 + this%w1))
            else if (z < this%z3) then
                grho_de = grho_de_today * fac2 * a**(-3 * (1 + this%w2))
            else
                grho_de = grho_de_today * fac3
            end if    
        else if (this%model == 4) then
            ! 5 Bins w
            fac1 = (1.0+this%z1)**(3 * (this%w0 - this%w1))
            fac2 = fac1 * (1.0+this%z2)**(3 * (this%w1 - this%w2))
            fac3 = fac2 * (1.0+this%z3)**(3 * (this%w2 - this%w3))
            fac4 = fac3 * (1.0+this%z4)**(3 * (this%w3 - this%w4))
            fac5 = fac4 * (1.0+this%z5)**(3 * (this%w4 - (-1)))

            if (z < this%z1) then
                grho_de = grho_de_today * a**(-3 * (1 + this%w0))
            else if (z < this%z2) then
                grho_de = grho_de_today * fac1 * a**(-3 * (1 + this%w1))
            else if (z < this%z3) then
                grho_de = grho_de_today * fac2 * a**(-3 * (1 + this%w2))
            else if (z < this%z4) then
                grho_de = grho_de_today * fac3 * a**(-3 * (1 + this%w3))
            else if (z < this%z5) then
                grho_de = grho_de_today * fac4 * a**(-3 * (1 + this%w4))
            else
                grho_de = grho_de_today * fac5
            end if    
        else if (this%model == 5) then
            ! 10 Bins w
            fac1 = (1.0+this%z1)**(3 * (this%w0 - this%w1))
            fac2 = fac1 * (1.0+this%z2)**(3 * (this%w1 - this%w2))
            fac3 = fac2 * (1.0+this%z3)**(3 * (this%w2 - this%w3))
            fac4 = fac3 * (1.0+this%z4)**(3 * (this%w3 - this%w4))
            fac5 = fac4 * (1.0+this%z5)**(3 * (this%w4 - this%w5))
            fac6 = fac5 * (1.0+this%z6)**(3 * (this%w5 - this%w6))
            fac7 = fac6 * (1.0+this%z7)**(3 * (this%w6 - this%w7))
            fac8 = fac7 * (1.0+this%z8)**(3 * (this%w7 - this%w8))
            fac9 = fac8 * (1.0+this%z9)**(3 * (this%w8 - this%w9))
            fac10 = fac9 * (1.0+this%z10)**(3 * (this%w9 - (-1)))

            if (z < this%z1) then
                grho_de = grho_de_today * a**(-3 * (1 + this%w0))
            else if (z < this%z2) then
                grho_de = grho_de_today * fac1 * a**(-3 * (1 + this%w1))
            else if (z < this%z3) then
                grho_de = grho_de_today * fac2 * a**(-3 * (1 + this%w2))
            else if (z < this%z4) then
                grho_de = grho_de_today * fac3 * a**(-3 * (1 + this%w3))
            else if (z < this%z5) then
                grho_de = grho_de_today * fac4 * a**(-3 * (1 + this%w4))
            else if (z < this%z6) then
                grho_de = grho_de_today * fac5 * a**(-3 * (1 + this%w5))
            else if (z < this%z7) then
                grho_de = grho_de_today * fac6 * a**(-3 * (1 + this%w6))
            else if (z < this%z8) then
                grho_de = grho_de_today * fac7 * a**(-3 * (1 + this%w7))
            else if (z < this%z9) then
                grho_de = grho_de_today * fac8 * a**(-3 * (1 + this%w8)) 
            else if (z < this%z10) then
                grho_de = grho_de_today * fac9 * a**(-3 * (1 + this%w9))
            else
                grho_de = grho_de_today * fac10
            end if    
        else if (this%model == 6) then            
            !2 Bins linear w(z)                       
            Delta_z1 = this%z1 - this%z0
            Delta_z2 = this%z2 - this%z1
            
            Delta_w1 = this%w1 - this%w0
            Delta_w2 = (-1.0)  - this%w1

            wa0 = Delta_w1/Delta_z1
            wa1 = Delta_w2/Delta_z2

            alpha0 = 3*(1.0+this%w0 - wa0*(1.0+this%z0))      
            alpha1 = 3*(1.0+this%w1 - wa1*(1.0+this%z1))
            
            z = 1._dl/a - 1
            if (z < this%z1) then
                grho_de = grho_de_today * ((1.0+z)/(1.0+this%z0))**alpha0*exp(3*wa0*(z-this%z0))
            else if (z < this%z2) then
                grho_de = grho_de_today * ((1.0+this%z1)/(1.0+this%z0))**alpha0*exp(3*wa0*(this%z1-this%z0)) * &
                                          ((1.0+z)/(1.0+this%z1))**alpha1*exp(3*wa1*(z -this%z1))
            else
                grho_de = grho_de_today * ((1.0+this%z1)/(1.0+this%z0))**alpha0*exp(3*wa0*(this%z1-this%z0)) * &
                                          ((1.0+this%z2)/(1.0+this%z1))**alpha1*exp(3*wa1*(this%z2-this%z1))
            end if 
        else if(this%model == 7) then 
            ! 2 bins quadratic w(z)            
            Delta_z1 = this%z1 - this%z0
            Delta_w1 = this%w1 - this%w0 
            Delta_z2 = this%z2 - this%z1
            !VM: LAST BIN IS ALWAYS -1.0
            Delta_w2 = (-1.0) - this%w1    

            wa0  = 2._dl*(Delta_w1/Delta_z1 - Delta_w2/Delta_z2)
            waa0 = -Delta_w1/Delta_z1**2 + 2._dl*Delta_w2/(Delta_z1*Delta_z2)
            wa1  =  2._dl*Delta_w2/Delta_z2
            waa1 = -Delta_w2/Delta_z2**2 
            wa2  = 0
            waa2 = 0

            A00 = 3._dl*(1.0+this%w0 -wa0*this%z0+waa0*this%z0**2._dl)
            A10 = 3._dl*(wa0-2*waa0*this%z0)
            A20 = 3._dl*waa0

            A01 = 3._dl*(1.0+this%w1-wa1*this%z1+waa1*this%z1**2._dl)
            A11 = 3._dl*(wa1-2*waa1*this%z1)
            A21 = 3*waa1
            
            if (z < this%z1) then
                grho_de = grho_de_today * &
                    (((1.0+ z)/(1.0+this%z0))**(A00-A10+A20))*exp((A10-A20)*(z-this%z0)+A20*(z**2-this%z0**2)/2)
            else if (z < this%z2) then
                grho_de = grho_de_today * &
                    (((1.0+this%z1)/(1.0+this%z0))**(A00-A10+A20))*exp((A10-A20)*(this%z1-this%z0)+A20*(this%z1**2-this%z0**2)/2) * &
                    (((1.0+z)/(1.0+this%z1))**(A01-A11+A21))*exp((A11-A21)*(z-this%z1)+A21*( z**2-this%z1**2)/2)
            else
                grho_de = grho_de_today * & 
                    (((1.0+this%z1)/(1.0+this%z0))**(A00-A10+A20))*exp((A10-A20)*(this%z1-this%z0)+A20*(this%z1**2-this%z0**2)/2) * &
                    (((1.0+this%z2)/(1.0+this%z1))**(A01-A11+A21))*exp((A11-A21)*(this%z2-this%z1)+A21*(this%z2**2-this%z1**2)/2)        
            end if            
        else 
            stop "[Late Fluid DE @TLateDE_grho_de] Invalid Dark Energy Model"
        end if
    end function TLateDE_grho_de

    subroutine TLateDE_Init(this, State)
        use classes
        use results
        class(TLateDE), intent(inout) :: this
        class(TCAMBdata), intent(in), target :: State

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
            stop "[Late Fluid DE @TLateDE_Effective_w_wa] Invalid Dark Energy Model (TLateDE_Effective_w_wa)"
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