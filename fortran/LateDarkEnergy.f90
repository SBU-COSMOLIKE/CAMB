module LateDE
    use DarkEnergyInterface
    use results
    use constants
    use classes
    implicit none

    private
    real(dl) :: grho_de_today, Integrate_tanh
    type, extends(TDarkEnergyModel) :: TLateDE
        integer  :: DEmodel
        real(dl) :: w0,w1,w2,w3,w4,w5,w6,w7,w8,w9
        real(dl) :: z1,z2,z3,z4,z5,z6,z7,z8,z9,z10
        real(dl) :: zc, sigma ! JVR: added tanh parameters
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
        real(dl) :: wa0,waa0,waaa0, wa1,waa1,waaa1, wa2,waa2,waaa2, wa3,waa3,waaa3, wa4,waa4,waaa4, a1, a2
        real(dl) :: Delta_z1, Delta_z2, Delta_z3, Delta_z4, Delta_z5
        real(dl) :: Delta_w1, Delta_w2, Delta_w3, Delta_w4, Delta_w5
        integer  :: i

        w_de = 0
        z = 1.0_dl/a - 1.0_dl

        if (this%DEmodel == 1) then
            ! Constant w
            w_de = this%w0
        else if (this%DEmodel == 2) then
            ! CPL parametrization w0wa
            w_de = this%w0 + this%w1*(1._dl - a)
        else if (this%DEmodel == 3) then
            ! Constant w: 2 bins
            if (z < this%z1) then
                w_de = this%w0
            else if (z < this%z2) then
                w_de = this%w1
            else
                w_de = -1.0_dl
            end if 
        else if (this%DEmodel == 4) then
            ! Constant w: 3 bins
            if (z < this%z1) then
                w_de = this%w0
            else if (z < this%z2) then
                w_de = this%w1
            else if (z < this%z3) then
                w_de = this%w2
            else
                w_de = -1.0_dl
            end if    
        else if (this%DEmodel == 5) then
            ! Constant w: 5 bins
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
                w_de = -1.0_dl
            end if
        else if (this%DEmodel == 6) then
            ! Linear w(z): 2 bins
            
            Delta_z1 = this%z1-0
            Delta_z2 = this%z2-this%z1
            Delta_w1 = this%w1-this%w0
            Delta_w2 = -1.0_dl-this%w1

            wa0 = Delta_w1/Delta_z1
            wa1 = Delta_w2/Delta_z2
            if (z < this%z1) then
                w_de = this%w0 + wa0 * (z - 0.0_dl)
            else if (z < this%z2) then
                w_de = this%w1 + wa1 * (z - this%z1)
            else
                w_de = -1.0_dl
            end if 
        else if (this%DEmodel == 7) then
            ! Linear w(z): 3 bins

            Delta_z1 = this%z1-0
            Delta_z2 = this%z2-this%z1
            Delta_z3 = this%z3-this%z2
            Delta_w1 = this%w1-this%w0
            Delta_w2 = this%w2-this%w1
            Delta_w3 = -1.0_dl-this%w2 

            wa0 = Delta_w1/Delta_z1
            wa1 = Delta_w2/Delta_z2
            wa2 = Delta_w3/Delta_z3
            if (z < this%z1) then
                w_de = this%w0 + wa0 * (z - 0.0_dl)
            else if (z < this%z2) then
                w_de = this%w1 + wa1 * (z - this%z1)
            else if (z < this%z3) then
                w_de = this%w2 + wa2 * (z - this%z2)
            else
                w_de = -1.0_dl    
            end if 
        else if (this%DEmodel == 8) then
            ! Linear w(z): 5 bins

            Delta_z1 = this%z1-0
            Delta_z2 = this%z2-this%z1
            Delta_z3 = this%z3-this%z2
            Delta_z4 = this%z4-this%z3
            Delta_z5 = this%z5-this%z4

            Delta_w1 = this%w1-this%w0
            Delta_w2 = this%w2-this%w1
            Delta_w3 = this%w3-this%w2
            Delta_w4 = this%w4-this%w3
            Delta_w5 = -1.0_dl-this%w4 

            wa0 = Delta_w1/Delta_z1
            wa1 = Delta_w2/Delta_z2
            wa2 = Delta_w3/Delta_z3
            wa3 = Delta_w4/Delta_z4
            wa4 = Delta_w5/Delta_z5
            if (z < this%z1) then
                w_de = this%w0 + wa0 * (z - 0.0_dl)
            else if (z < this%z2) then
                w_de = this%w1 + wa1 * (z - this%z1)
            else if (z < this%z3) then
                w_de = this%w2 + wa2 * (z - this%z2)
            else if (z < this%z4) then
                w_de = this%w3 + wa3 * (z - this%z3)
            else if (z < this%z5) then
                w_de = this%w4 + wa4 * (z - this%z4)
            else
                w_de = -1.0_dl    
            end if                          
        else if (this%DEmodel == 9) then
            ! Quadratic w(z): 2 bins
            ! Boundary conditions -see my notes-

            Delta_z1 = this%z1-0
            Delta_z2 = this%z2-this%z1
            Delta_w1 = this%w1-this%w0
            Delta_w2 = -1.0_dl-this%w1

            wa0  = 2._dl*(Delta_w1/Delta_z1 - Delta_w2/Delta_z2)
            waa0 = -Delta_w1/Delta_z1**2 + 2._dl*Delta_w2/(Delta_z1*Delta_z2)
            wa1  =  2._dl*Delta_w2/Delta_z2
            waa1 =  -Delta_w2/Delta_z2**2
            ! Equation of state
            if (z < this%z1) then
                w_de = this%w0 + wa0 * (z-0) + waa0*(z-0.0_dl)**2
            else if (z < this%z2) then
                w_de = this%w1 + wa1 * (z-this%z1) + waa1*(z-this%z1)**2
            else
                w_de = -1.0_dl
            end if
        else if (this%DEmodel == 10) then
            ! Quadratic w(z): 3 bins
            ! Boundary conditions -see my notes-

            Delta_z1 = this%z1-0
            Delta_z2 = this%z2-this%z1
            Delta_z3 = this%z3-this%z2
            Delta_w1 = this%w1-this%w0
            Delta_w2 = this%w2-this%w1
            Delta_w3 = -1.0_dl-this%w2 

            !See Mathematica Notebook - CAMBLateDE_boundary_conditions.nb
            wa0  = 2._dl*(Delta_w1/Delta_z1 - Delta_w2/Delta_z2 + Delta_w3/Delta_z3)
            waa0 = -Delta_w1/Delta_z1**2 + (2._dl/Delta_z1)*(Delta_w2/Delta_z2 - Delta_w3/Delta_z3)
            wa1  =  2._dl*(Delta_w2/Delta_z2 - Delta_w3/Delta_z3)
            waa1 =  -Delta_w2/Delta_z2**2 + 2._dl*Delta_w3/(Delta_z2*Delta_z3)
            wa2  = 2._dl*Delta_w3/Delta_z3
            waa2 = -Delta_w3/Delta_z3**2
            ! Equation of state
            if (z < this%z1) then
                w_de = this%w0 + wa0*(z-0.0_dl) + waa0*(z-0.0_dl)**2
            else if (z < this%z2) then
                w_de = this%w1 + wa1*(z-this%z1) + waa1*(z-this%z1)**2
            else if (z < this%z3) then
                w_de = this%w2 + wa2*(z-this%z2) + waa2*(z-this%z2)**2
            else
                w_de = -1.0_dl
            end if
        else if (this%DEmodel == 11) then
            ! Quadratic w(z): 5 bins
            ! Boundary conditions -see my notes-

            Delta_z1 = this%z1-0
            Delta_z2 = this%z2-this%z1
            Delta_z3 = this%z3-this%z2
            Delta_z4 = this%z4-this%z3
            Delta_z5 = this%z5-this%z4

            Delta_w1 = this%w1-this%w0
            Delta_w2 = this%w2-this%w1
            Delta_w3 = this%w3-this%w2
            Delta_w4 = this%w4-this%w3
            Delta_w5 = -1.0_dl-this%w4 

            !See Mathematica Notebook - CAMBLateDE_boundary_conditions.nb
            wa0  = 2._dl*(Delta_w1/Delta_z1 - Delta_w2/Delta_z2 + Delta_w3/Delta_z3 - Delta_w4/Delta_z4 + Delta_w5/Delta_z5)
            waa0 = -Delta_w1/Delta_z1**2 + (2._dl/Delta_z1)*(Delta_w2/Delta_z2 - Delta_w3/Delta_z3 + Delta_w4/Delta_z4 - Delta_w5/Delta_z5)
            
            wa1  = 2._dl*(Delta_w2/Delta_z2 - Delta_w3/Delta_z3 + Delta_w4/Delta_z4 - Delta_w5/Delta_z5)
            waa1 = -Delta_w2/Delta_z2**2 + (2._dl/Delta_z2)*(Delta_w3/Delta_z3 - Delta_w4/Delta_z4 + Delta_w5/Delta_z5)
            
            wa2  = 2._dl*(Delta_w3/Delta_z3 - Delta_w4/Delta_z4 + Delta_w5/Delta_z5)
            waa2 = -Delta_w3/Delta_z3**2 + (2._dl/Delta_z3)*(Delta_w4/Delta_z4 - Delta_w5/Delta_z5)

            wa3  = 2._dl*(Delta_w4/Delta_z4 - Delta_w5/Delta_z5)
            waa3 = -Delta_w4/Delta_z4**2 + (2._dl/Delta_z4)*(Delta_w5/Delta_z5)

            wa4  = 2._dl*(Delta_w5/Delta_z5)
            waa4 = -Delta_w5/Delta_z5**2

            ! Equation of state
            if (z < this%z1) then
                w_de = this%w0 + wa0*(z-0.0_dl) + waa0*(z-0.0_dl)**2
            else if (z < this%z2) then
                w_de = this%w1 + wa1*(z-this%z1) + waa1*(z-this%z1)**2
            else if (z < this%z3) then
                w_de = this%w2 + wa2*(z-this%z2) + waa2*(z-this%z2)**2
            else if (z < this%z4) then
                w_de = this%w3 + wa3*(z-this%z3) + waa3*(z-this%z3)**2
            else if (z < this%z5) then
                w_de = this%w4 + wa4*(z-this%z4) + waa4*(z-this%z4)**2
            else
                w_de = -1.0_dl
            end if             
        else if (this%DEmodel == 12) then
            ! Cubic w(z): 2 bins
            ! Boundary conditions -see my notes-

            Delta_z1 = this%z1-0
            Delta_z2 = this%z2-this%z1
            Delta_w1 = this%w1-this%w0
            Delta_w2 = -1.0_dl-this%w1

            wa0   = 3*Delta_w1/Delta_z1 - 3*Delta_w2/Delta_z2 * (2+Delta_z1/Delta_z2)
            waa0  = -3*Delta_w1/Delta_z1**2 + 3*Delta_w2/Delta_z2*(3/Delta_z1 + 2/Delta_z2)
            waaa0 = Delta_w1/Delta_z1**3 - 3*Delta_w2/(Delta_z1*Delta_z2) * (1/Delta_z1 + 1/Delta_z2)
            wa1   = 3*Delta_w2/Delta_z2
            waa1  = -3*Delta_w2/Delta_z2**2
            waaa1 = Delta_w2/Delta_z2**3
            ! Equation of state
            if (z < this%z1) then
                w_de = this%w0 + wa0 * (z-0) + waa0*(z-0)**2 + waaa0*(z-0)**3
            else if (z < this%z2) then
                w_de = this%w1 + wa1 * (z-this%z1) + waa1*(z-this%z1)**2 + waaa1*(z-this%z1)**3
            else
                w_de = -1.0_dl
            end if
        else if (this%DEmodel == 13) then
            ! Cubic w(z): 3 bins
            ! Boundary conditions -see my notes-

            Delta_z1 = this%z1-0
            Delta_z2 = this%z2-this%z1
            Delta_z3 = this%z3-this%z2
            Delta_w1 = this%w1-this%w0
            Delta_w2 = this%w2-this%w1
            Delta_w3 = -1.0_dl-this%w2 

            wa0   = 3*Delta_w1/Delta_z1 - 3*(Delta_w2/Delta_z2)*(2 + Delta_z1/Delta_z2) + 3*(Delta_w3/Delta_z3)*(4 + 3*Delta_z1/Delta_z2 + 2*Delta_z1/Delta_z3 + 2*Delta_z2/Delta_z3)
            waa0  = -3*Delta_w1/Delta_z1**2 + 3*(Delta_w2/Delta_z2)*(3/Delta_z1 + 2/Delta_z2) -3*(Delta_w3/Delta_z3)*(6/Delta_z1 + 6/Delta_z2 + 4/Delta_z3 + 3*Delta_z2/(Delta_z1*Delta_z3))
            waaa0 = Delta_w1/Delta_z3**3 - 3*Delta_w2/(Delta_z1*Delta_z2)*(1/Delta_z1 + 1/Delta_z2) + 3*Delta_w3/(Delta_z1*Delta_z3)*(2/Delta_z1 + 3/Delta_z2 + 2/Delta_z3 + Delta_z2/(Delta_z1*Delta_z3))
            wa1   = 3*Delta_w2/Delta_z2 - 3*Delta_w3/Delta_z3*(2 + Delta_z2/Delta_z3)
            waa1  = -3*Delta_w2/Delta_z2**2 + 3*Delta_w3/Delta_z3*(3/Delta_z2 + 2/Delta_z3)
            waaa1 = Delta_w2/Delta_z2**3 - 3*Delta_w3/(Delta_z2*Delta_z3)*(1/Delta_z2 + 1/Delta_z3)
            wa2   = 3*Delta_w3/Delta_z3
            waa2  = -3*Delta_w3/Delta_z3**2
            waaa2 = Delta_w3/Delta_z3**3
            ! Equation of state
            if (z < this%z1) then
                w_de = this%w0 + wa0 * (z-0) + waa0*(z-0)**2 + waaa0*(z-0)**3
            else if (z < this%z2) then
                w_de = this%w1 + wa1 * (z-this%z1) + waa1*(z-this%z1)**2 + waaa1*(z-this%z1)**3
            else if (z < this%z3) then
                w_de = this%w2 + wa2 * (z-this%z2) + waa2*(z-this%z2)**2 + waaa2*(z-this%z2)**3
            else 
                w_de = -1.0_dl    
            end if

        else if (this%DEmodel == 14) then
            ! Cubic w(z): 5 bins
            ! Boundary conditions -see my notes-

            Delta_z1 = this%z1-0
            Delta_z2 = this%z2-this%z1
            Delta_z3 = this%z3-this%z2
            Delta_z4 = this%z4-this%z3
            Delta_z5 = this%z5-this%z4

            Delta_w1 = this%w1-this%w0
            Delta_w2 = this%w2-this%w1
            Delta_w3 = this%w3-this%w2
            Delta_w4 = this%w4-this%w3
            Delta_w5 = -1.0_dl-this%w4 

            !See Mathematica Notebook - CAMBLateDE_boundary_conditions.nb
            wa0   = (3*(Delta_w1))/(Delta_z1) - (3*(Delta_w2)*(Delta_z1))/(Delta_z2)**2 - (6*(Delta_w2))/(Delta_z2) + &
                    ((Delta_w3)*(6*(Delta_z1) + 6*(Delta_z2)))/(Delta_z3)**2 + ((Delta_w3)*(12 + (9*(Delta_z1))/(Delta_z2)))/(Delta_z3) + &
                    ((Delta_w4)*(-12*(Delta_z1) - 12*(Delta_z2) - 12*(Delta_z3) - (9*(Delta_z1)*(Delta_z3))/(Delta_z2)))/(Delta_z4)**2 + &
                    ((Delta_w4)*(-24 - (18*(Delta_z1))/(Delta_z2) - (18*(Delta_z1))/(Delta_z3) - (18*(Delta_z2))/(Delta_z3)))/(Delta_z4) + & 
                    ((Delta_w5)*(24*(Delta_z1) + 24*(Delta_z2) + 24*(Delta_z3) + (18*(Delta_z1)*(Delta_z3))/(Delta_z2) + 24*(Delta_z4) + &
                    (18*(Delta_z1)*(Delta_z4))/(Delta_z2) + (18*(Delta_z1)*(Delta_z4))/(Delta_z3) + (18*(Delta_z2)*(Delta_z4))/(Delta_z3)))/(Delta_z5)**2 + &
                    ((Delta_w5)*(48 + (36*(Delta_z1))/(Delta_z2) + (36*(Delta_z1))/(Delta_z3) + (36*(Delta_z2))/(Delta_z3) + (36*(Delta_z1))/(Delta_z4) + &
                    (36*(Delta_z2))/(Delta_z4) + (36*(Delta_z3))/(Delta_z4) + (27*(Delta_z1)*(Delta_z3))/((Delta_z2)*(Delta_z4))))/(Delta_z5)
            waa0  = (-3*(Delta_w1))/(Delta_z1)**2 + (6*(Delta_w2))/(Delta_z2)**2 + (9*(Delta_w2))/((Delta_z1)*(Delta_z2)) + &
                    ((Delta_w3)*(-12 - (9*(Delta_z2))/(Delta_z1)))/(Delta_z3)**2 + ((Delta_w3)*(-18/(Delta_z1) - 18/(Delta_z2)))/(Delta_z3) + &
                    ((Delta_w4)*(24 + (18*(Delta_z2))/(Delta_z1) + (18*(Delta_z3))/(Delta_z1) + (18*(Delta_z3))/(Delta_z2)))/(Delta_z4)**2 + &
                    ((Delta_w4)*(36/(Delta_z1) + 36/(Delta_z2) + 36/(Delta_z3) + (27*(Delta_z2))/((Delta_z1)*(Delta_z3))))/(Delta_z4) + &
                    ((Delta_w5)*(-48 - (36*(Delta_z2))/(Delta_z1) - (36*(Delta_z3))/(Delta_z1) - (36*(Delta_z3))/(Delta_z2) - (36*(Delta_z4))/(Delta_z1) - &
                    (36*(Delta_z4))/(Delta_z2) - (36*(Delta_z4))/(Delta_z3) - (27*(Delta_z2)*(Delta_z4))/((Delta_z1)*(Delta_z3))))/(Delta_z5)**2 + &
                    ((Delta_w5)*(-72/(Delta_z1) - 72/(Delta_z2) - 72/(Delta_z3) - (54*(Delta_z2))/((Delta_z1)*(Delta_z3)) - 72/(Delta_z4) - &
                    (54*(Delta_z2))/((Delta_z1)*(Delta_z4)) - (54*(Delta_z3))/((Delta_z1)*(Delta_z4)) - (54*(Delta_z3))/((Delta_z2)*(Delta_z4))))/(Delta_z5)
            waaa0 = (Delta_w1)/(Delta_z1)**3 - (3*(Delta_w2))/((Delta_z1)*(Delta_z2)**2) - (3*(Delta_w2))/((Delta_z1)**2*(Delta_z2)) + &
                    ((Delta_w3)*(6/(Delta_z1) + (3*(Delta_z2))/(Delta_z1)**2))/(Delta_z3)**2 + ((Delta_w3)*(6/(Delta_z1)**2 + 9/((Delta_z1)*(Delta_z2))))/(Delta_z3) + &
                    ((Delta_w4)*(-12/(Delta_z1) - (6*(Delta_z2))/(Delta_z1)**2 - (6*(Delta_z3))/(Delta_z1)**2 - (9*(Delta_z3))/((Delta_z1)*(Delta_z2))))/&
                    (Delta_z4)**2 + ((Delta_w4)*(-12/(Delta_z1)**2 - 18/((Delta_z1)*(Delta_z2)) - 18/((Delta_z1)*(Delta_z3)) - &
                    (9*(Delta_z2))/((Delta_z1)**2*(Delta_z3))))/(Delta_z4) + &
                    ((Delta_w5)*(24/(Delta_z1) + (12*(Delta_z2))/(Delta_z1)**2 + (12*(Delta_z3))/(Delta_z1)**2 + (18*(Delta_z3))/((Delta_z1)*(Delta_z2)) + &
                    (12*(Delta_z4))/(Delta_z1)**2 + (18*(Delta_z4))/((Delta_z1)*(Delta_z2)) + (18*(Delta_z4))/((Delta_z1)*(Delta_z3)) + &
                    (9*(Delta_z2)*(Delta_z4))/((Delta_z1)**2*(Delta_z3))))/(Delta_z5)**2 + &
                    ((Delta_w5)*(24/(Delta_z1)**2 + 36/((Delta_z1)*(Delta_z2)) + 36/((Delta_z1)*(Delta_z3)) + (18*(Delta_z2))/((Delta_z1)**2*(Delta_z3)) + &
                    36/((Delta_z1)*(Delta_z4)) + (18*(Delta_z2))/((Delta_z1)**2*(Delta_z4)) + (18*(Delta_z3))/((Delta_z1)**2*(Delta_z4)) + &
                    (27*(Delta_z3))/((Delta_z1)*(Delta_z2)*(Delta_z4))))/(Delta_z5)
            
            wa1   = (3*(Delta_w2))/(Delta_z2) - (3*(Delta_w3)*(Delta_z2))/(Delta_z3)**2 - (6*(Delta_w3))/(Delta_z3) + &
                    ((Delta_w4)*(6*(Delta_z2) + 6*(Delta_z3)))/(Delta_z4)**2 + ((Delta_w4)*(12 + (9*(Delta_z2))/(Delta_z3)))/(Delta_z4) + &
                    ((Delta_w5)*(-12*(Delta_z2) - 12*(Delta_z3) - 12*(Delta_z4) - (9*(Delta_z2)*(Delta_z4))/(Delta_z3)))/(Delta_z5)**2 + &
                    ((Delta_w5)*(-24 - (18*(Delta_z2))/(Delta_z3) - (18*(Delta_z2))/(Delta_z4) - (18*(Delta_z3))/(Delta_z4)))/(Delta_z5)
            waa1  = (-3*(Delta_w2))/(Delta_z2)**2 + (6*(Delta_w3))/(Delta_z3)**2 + (9*(Delta_w3))/((Delta_z2)*(Delta_z3)) + &
                    ((Delta_w4)*(-12 - (9*(Delta_z3))/(Delta_z2)))/(Delta_z4)**2 + ((Delta_w4)*(-18/(Delta_z2) - 18/(Delta_z3)))/(Delta_z4) + &
                    ((Delta_w5)*(24 + (18*(Delta_z3))/(Delta_z2) + (18*(Delta_z4))/(Delta_z2) + (18*(Delta_z4))/(Delta_z3)))/(Delta_z5)**2 + &
                    ((Delta_w5)*(36/(Delta_z2) + 36/(Delta_z3) + 36/(Delta_z4) + (27*(Delta_z3))/((Delta_z2)*(Delta_z4))))/(Delta_z5)
            waaa1 = (Delta_w2)/(Delta_z2)**3 - (3*(Delta_w3))/((Delta_z2)*(Delta_z3)**2) - (3*(Delta_w3))/((Delta_z2)**2*(Delta_z3)) + &
                    ((Delta_w4)*(6/(Delta_z2) + (3*(Delta_z3))/(Delta_z2)**2))/(Delta_z4)**2 + ((Delta_w4)*(6/(Delta_z2)**2 + 9/((Delta_z2)*(Delta_z3))))/(Delta_z4) + &
                    ((Delta_w5)*(-12/(Delta_z2) - (6*(Delta_z3))/(Delta_z2)**2 - (6*(Delta_z4))/(Delta_z2)**2 - (9*(Delta_z4))/((Delta_z2)*(Delta_z3))))/&
                    (Delta_z5)**2 + ((Delta_w5)*(-12/(Delta_z2)**2 - 18/((Delta_z2)*(Delta_z3)) - 18/((Delta_z2)*(Delta_z4)) - &
                    (9*(Delta_z3))/((Delta_z2)**2*(Delta_z4))))/(Delta_z5)
            
            wa2   = (3*(Delta_w3))/(Delta_z3) - (3*(Delta_w4)*(Delta_z3))/(Delta_z4)**2 - (6*(Delta_w4))/(Delta_z4) + &
                    ((Delta_w5)*(6*(Delta_z3) + 6*(Delta_z4)))/(Delta_z5)**2 + ((Delta_w5)*(12 + (9*(Delta_z3))/(Delta_z4)))/(Delta_z5)
            waa2  = (-3*(Delta_w3))/(Delta_z3)**2 + (6*(Delta_w4))/(Delta_z4)**2 + (9*(Delta_w4))/((Delta_z3)*(Delta_z4)) + &
                    ((Delta_w5)*(-12 - (9*(Delta_z4))/(Delta_z3)))/(Delta_z5)**2 + ((Delta_w5)*(-18/(Delta_z3) - 18/(Delta_z4)))/(Delta_z5)
            waaa2 = (Delta_w3)/(Delta_z3)**3 - (3*(Delta_w4))/((Delta_z3)*(Delta_z4)**2) - (3*(Delta_w4))/((Delta_z3)**2*(Delta_z4)) + &
                    ((Delta_w5)*(6/(Delta_z3) + (3*(Delta_z4))/(Delta_z3)**2))/(Delta_z5)**2 + ((Delta_w5)*(6/(Delta_z3)**2 + 9/((Delta_z3)*(Delta_z4))))/(Delta_z5)

            wa3   = (3*(Delta_w4))/(Delta_z4) - (3*(Delta_w5)*(Delta_z4))/(Delta_z5)**2 - (6*(Delta_w5))/(Delta_z5)
            waa3  = (-3*(Delta_w4))/(Delta_z4)**2 + (6*(Delta_w5))/(Delta_z5)**2 + (9*(Delta_w5))/((Delta_z4)*(Delta_z5))
            waaa3 = (Delta_w4)/(Delta_z4)**3 - (3*(Delta_w5))/((Delta_z4)*(Delta_z5)**2) - (3*(Delta_w5))/((Delta_z4)**2*(Delta_z5))

            wa4   = (3*(Delta_w5))/(Delta_z5)
            waa4  = (-3*(Delta_w5))/(Delta_z5)**2
            waaa4 = (Delta_w5)/(Delta_z5)**3

            ! Equation of state
            if (z < this%z1) then
                w_de = this%w0 + wa0 * (z-0) + waa0*(z-0)**2 + waaa0*(z-0)**3
            else if (z < this%z2) then
                w_de = this%w1 + wa1 * (z-this%z1) + waa1*(z-this%z1)**2 + waaa1*(z-this%z1)**3
            else if (z < this%z3) then
                w_de = this%w2 + wa2 * (z-this%z2) + waa2*(z-this%z2)**2 + waaa2*(z-this%z2)**3
            else if (z < this%z4) then
                w_de = this%w3 + wa3 * (z-this%z3) + waa3*(z-this%z3)**2 + waaa3*(z-this%z3)**3
            else if (z < this%z5) then    
                w_de = this%w4 + wa4 * (z-this%z4) + waa4*(z-this%z4)**2 + waaa4*(z-this%z4)**3
            else 
                w_de = -1.0_dl
            end if            
        else if (this%DEmodel == 15) then
            ! tanh model - w(z) transition
            if (z < this%z1 - this%sigma) then
                w_de = this%w0
            else if (z > this%z1 + this%sigma) then
                w_de = this%w1
            else
                a1 = (this%w1 - this%w0)/this%sigma/2._dl
                a2 = this%w1 - a1*(this%z1 + this%sigma)
                w_de = a2 + a1*z
            end if
            !DHFS MOD TANH START
        else if (this%DEmodel == 16) then 
                ! Numeric tanh
                w_de = this%w0 + (this%w1-this%w0)/2._dl * ( 1.0_dl + tanh((z - this%z1)/this%sigma) )
            !DHFS MOD TANH END
        else        
            stop "[Late Fluid DE @TLateDE_w_de] Invalid Dark Energy Model"   
        end if
    end function TLateDE_w_de

    !DHFS MOD TANH START
    function kernel_tanh(this,z)
        class(TLateDE) :: this
        real(dl) :: kernel_tanh, w_de
        real(dl), intent(in) :: z

        w_de = this%w0 + (this%w1-this%w0)/2._dl * ( 1.0_dl + tanh((z-this%z1)/this%sigma) )
        kernel_tanh = (1.0_dl + w_de) / (1.0_dl+z)
    end function kernel_tanh
    !DHFS MOD TANH END

    function TLateDE_grho_de(this, a) result(grho_de)
        ! Returns 8*pi*G * rho_de, no factor of a^4
        class(TLateDE) :: this
        real(dl), intent(in) :: a
        real(dl) :: grho_de, z, a1, a2, za, zb
        real(dl) :: alpha0,alpha1,alpha2,alpha3,alpha4
        real(dl) :: Delta_z1, Delta_z2, Delta_z3, Delta_z4, Delta_z5
        real(dl) :: Delta_w1, Delta_w2, Delta_w3, Delta_w4, Delta_w5  
        real(dl) :: fac1,fac2,fac3,fac4,fac5,fac6,fac7,fac8,fac9,fac10
        real(dl) :: wa0,waa0,waaa0, A00,A10,A20,A30 ! factors for the 1st bin
        real(dl) :: wa1,waa1,waaa1, A01,A11,A21,A31 ! factors for the 2st bin
        real(dl) :: wa2,waa2,waaa2, A02,A12,A22,A32 ! factors for the 3st bin
        real(dl) :: wa3,waa3,waaa3, A03,A13,A23,A33 ! factors for the 4th bin
        real(dl) :: wa4,waa4,waaa4, A04,A14,A24,A34 ! factors for the 5th bin
        integer :: i
        integer, parameter :: npts = 50

        grho_de = 0
        z = 1.0_dl/a - 1.0_dl

        if (this%DEmodel == 1) then
            ! w constant
            grho_de = grho_de_today * a**(-3 * (1 + this%w0))
        else if (this%DEmodel == 2) then
            ! CPL w0-wa
            grho_de = grho_de_today * a**(-3 * (1 + this%w0 + this%w1)) * exp(-3 * this%w1 * (1 - a))
        else if (this%DEmodel == 3) then
            ! Constant w: 2 bins
            fac1 = (1.0_dl+this%z1)**(3.0_dl * (this%w0 - this%w1))
            fac2 = fac1 * (1.0_dl+this%z2)**(3.0_dl * (this%w1 - (-1.0_dl)))      
            if (z < this%z1) then
                grho_de = grho_de_today * a**(-3.0_dl * (1.0_dl + this%w0))
            else if (z < this%z2) then
                grho_de = grho_de_today * fac1 * a**(-3.0_dl * (1.0_dl + this%w1))
            else
                grho_de = grho_de_today * fac2
            end if  
        else if (this%DEmodel == 4) then
            ! Constant w: 3 bins
            fac1 = (1.0_dl+this%z1)**(3.0_dl * (this%w0 - this%w1))
            fac2 = fac1 * (1.0_dl+this%z2)**(3.0_dl * (this%w1 - this%w2))
            fac3 = fac2 * (1.0_dl+this%z3)**(3.0_dl * (this%w2 - (-1.0_dl)))      
            if (z < this%z1) then
                grho_de = grho_de_today * a**(-3.0_dl * (1.0_dl + this%w0))
            else if (z < this%z2) then
                grho_de = grho_de_today * fac1 * a**(-3.0_dl * (1.0_dl + this%w1))
            else if (z < this%z3) then
                grho_de = grho_de_today * fac2 * a**(-3.0_dl * (1.0_dl + this%w2))
            else
                grho_de = grho_de_today * fac3
            end if    
        else if (this%DEmodel == 5) then
            ! Constant w: 5 bins
            fac1 = (1.0_dl+this%z1)**(3.0_dl * (this%w0 - this%w1))
            fac2 = fac1 * (1.0_dl+this%z2)**(3.0_dl * (this%w1 - this%w2))
            fac3 = fac2 * (1.0_dl+this%z3)**(3.0_dl * (this%w2 - this%w3))
            fac4 = fac3 * (1.0_dl+this%z4)**(3.0_dl * (this%w3 - this%w4))
            fac5 = fac4 * (1.0_dl+this%z5)**(3.0_dl * (this%w4 - (-1.0_dl)))            
            if (z < this%z1) then
                grho_de = grho_de_today * a**(-3.0_dl * (1.0_dl + this%w0))
            else if (z < this%z2) then
                grho_de = grho_de_today * fac1 * a**(-3.0_dl * (1.0_dl + this%w1))
            else if (z < this%z3) then
                grho_de = grho_de_today * fac2 * a**(-3.0_dl * (1.0_dl + this%w2))
            else if (z < this%z4) then
                grho_de = grho_de_today * fac3 * a**(-3.0_dl * (1.0_dl + this%w3))
            else if (z < this%z5) then
                grho_de = grho_de_today * fac4 * a**(-3.0_dl * (1.0_dl + this%w4))
            else
                grho_de = grho_de_today * fac5
            end if    
        else if (this%DEmodel == 6) then
            ! Linear w(z): 2 bins

            Delta_z1 = this%z1-0
            Delta_z2 = this%z2-this%z1
            Delta_w1 = this%w1-this%w0
            Delta_w2 = -1.0_dl-this%w1

            wa0 = Delta_w1/Delta_z1
            wa1 = Delta_w2/Delta_z2

            alpha0 = 3.0_dl*(1.0_dl+this%w0-wa0*(1.0_dl+0))      
            alpha1 = 3.0_dl*(1.0_dl+this%w1-wa1*(1.0_dl+this%z1))

            fac1 = ((1.0_dl+this%z1)/(1.0_dl+0))**alpha0*exp(3.0_dl*wa0*(this%z1-0))
            fac2 = ((1.0_dl+this%z2)/(1.0_dl+this%z1))**alpha1*exp(3.0_dl*wa1*(this%z2-this%z1))

            if (z < this%z1) then
                grho_de = grho_de_today * ((1.0_dl+z )/(1.0_dl+0))**alpha0*exp(3.0_dl*wa0*(z-0))
            else if (z < this%z2) then
                grho_de = grho_de_today * fac1 * ((1.0_dl+z)/(1.0_dl+this%z1))**alpha1*exp(3.0_dl*wa1*(z -this%z1))
            else
                grho_de = grho_de_today * fac1 * fac2
            end if 
        else if (this%DEmodel == 7) then
            ! Linear w(z): 3 bins

            Delta_z1 = this%z1-0
            Delta_z2 = this%z2-this%z1
            Delta_z3 = this%z3-this%z2
            Delta_w1 = this%w1-this%w0
            Delta_w2 = this%w2-this%w1
            Delta_w3 = -1.0_dl-this%w2 

            wa0 = Delta_w1/Delta_z1
            wa1 = Delta_w2/Delta_z2
            wa2 = Delta_w3/Delta_z3

            alpha0 = 3.0_dl*(1.0_dl+this%w0-wa0*(1.0_dl+0))      
            alpha1 = 3.0_dl*(1.0_dl+this%w1-wa1*(1.0_dl+this%z1))
            alpha2 = 3.0_dl*(1.0_dl+this%w2-wa2*(1.0_dl+this%z2))

            fac1 = ((1.0_dl+this%z1)/(1.0_dl+0))**alpha0*exp(3.0_dl*wa0*(this%z1-0))
            fac2 = ((1.0_dl+this%z2)/(1.0_dl+this%z1))**alpha1*exp(3.0_dl*wa1*(this%z2-this%z1))
            fac3 = ((1.0_dl+this%z3)/(1.0_dl+this%z2))**alpha2*exp(3.0_dl*wa2*(this%z3-this%z2))

            if (z < this%z1) then
                grho_de = grho_de_today * ((1.0_dl+z )/(1.0_dl+0))**alpha0*exp(3.0_dl*wa0*(z-0))
            else if (z < this%z2) then
                grho_de = grho_de_today * fac1 * ((1.0_dl+z)/(1.0_dl+this%z1))**alpha1*exp(3.0_dl*wa1*(z-this%z1))
            else if (z < this%z3) then
                grho_de = grho_de_today * fac1 * fac2 * ((1.0_dl+z)/(1.0_dl+this%z2))**alpha2*exp(3.0_dl*wa2*(z-this%z2))
            else
                grho_de = grho_de_today * fac1 * fac2 * fac3
                                        
            end if 
        else if (this%DEmodel == 8) then
            ! Linear w(z): 5 bins

            Delta_z1 = this%z1-0
            Delta_z2 = this%z2-this%z1
            Delta_z3 = this%z3-this%z2
            Delta_z4 = this%z4-this%z3
            Delta_z5 = this%z5-this%z4

            Delta_w1 = this%w1-this%w0
            Delta_w2 = this%w2-this%w1
            Delta_w3 = this%w3-this%w2
            Delta_w4 = this%w4-this%w3
            Delta_w5 = -1.0_dl-this%w4 

            wa0 = Delta_w1/Delta_z1
            wa1 = Delta_w2/Delta_z2
            wa2 = Delta_w3/Delta_z3
            wa3 = Delta_w4/Delta_z4
            wa4 = Delta_w5/Delta_z5

            alpha0 = 3.0_dl*(1.0_dl+this%w0-wa0*(1.0_dl+0))      
            alpha1 = 3.0_dl*(1.0_dl+this%w1-wa1*(1.0_dl+this%z1))
            alpha2 = 3.0_dl*(1.0_dl+this%w2-wa2*(1.0_dl+this%z2))
            alpha3 = 3.0_dl*(1.0_dl+this%w3-wa3*(1.0_dl+this%z3))
            alpha4 = 3.0_dl*(1.0_dl+this%w3-wa4*(1.0_dl+this%z4))

            fac1 = ((1.0_dl+this%z1)/(1.0_dl+0))**alpha0*exp(3.0_dl*wa0*(this%z1-0))
            fac2 = ((1.0_dl+this%z2)/(1.0_dl+this%z1))**alpha1*exp(3.0_dl*wa1*(this%z2-this%z1))
            fac3 = ((1.0_dl+this%z3)/(1.0_dl+this%z2))**alpha2*exp(3.0_dl*wa2*(this%z3-this%z2))
            fac4 = ((1.0_dl+this%z4)/(1.0_dl+this%z3))**alpha3*exp(3.0_dl*wa3*(this%z4-this%z3))
            fac5 = ((1.0_dl+this%z5)/(1.0_dl+this%z4))**alpha4*exp(3.0_dl*wa4*(this%z5-this%z4))

            if (z < this%z1) then
                grho_de = grho_de_today * ((1.0_dl+z )/(1.0_dl+0))**alpha0*exp(3.0_dl*wa0*(z-0))
            else if (z < this%z2) then
                grho_de = grho_de_today * fac1 * ((1.0_dl+z)/(1.0_dl+this%z1))**alpha1*exp(3.0_dl*wa1*(z-this%z1))
            else if (z < this%z3) then
                grho_de = grho_de_today * fac1 * fac2 * ((1.0_dl+z)/(1.0_dl+this%z2))**alpha2*exp(3.0_dl*wa2*(z-this%z2))
            else if (z < this%z4) then
                grho_de = grho_de_today * fac1 * fac2 * fac3 * ((1.0_dl+z)/(1.0_dl+this%z3))**alpha3*exp(3.0_dl*wa3*(z-this%z3))
            else if (z < this%z5) then
                grho_de = grho_de_today * fac1 * fac2 * fac3 * fac4 * ((1.0_dl+z)/(1.0_dl+this%z4))**alpha4*exp(3.0_dl*wa4*(z-this%z4))
            else
                grho_de = grho_de_today * fac1 * fac2 * fac3 * fac4 * fac5
            end if                         
        else if(this%DEmodel == 9) then 
            ! Quadratic w(z): 2 bins  

            Delta_z1 = this%z1-0
            Delta_z2 = this%z2-this%z1
            Delta_w1 = this%w1-this%w0
            Delta_w2 = -1.0_dl-this%w1

            wa0  = 2._dl*(Delta_w1/Delta_z1 - Delta_w2/Delta_z2)
            waa0 = -Delta_w1/Delta_z1**2 + 2._dl*Delta_w2/(Delta_z1*Delta_z2)
            wa1  = 2._dl*Delta_w2/Delta_z2
            waa1 = -Delta_w2/Delta_z2**2 

            A00 = 3._dl*(1._dl + this%w0 - wa0*0 + waa0*0**2._dl)
            A10 = 3._dl*(wa0 - 2*waa0*0)
            A20 = 3._dl*waa0

            A01 = 3._dl*(1 + this%w1 - wa1*this%z1 + waa1*this%z1**2._dl)
            A11 = 3._dl*(wa1 - 2*waa1*this%z1)
            A21 = 3._dl*waa1

            fac1 = (((1+this%z1)/(1+0))**(A00-A10+A20))*exp((A10-A20)*(this%z1-0)+A20*(this%z1**2-0**2)/2)
            fac2 = (((1+this%z2)/(1+this%z1))**(A01-A11+A21))*exp((A11-A21)*(this%z2-this%z1)+A21*(this%z2**2-this%z1**2)/2)

            if (z < this%z1) then
                grho_de = grho_de_today * (((1+z)/(1+0))**(A00-A10+A20))*exp((A10-A20)*(z-0)+A20*(z**2-0**2)/2)
            else if (z < this%z2) then
                grho_de = grho_de_today * fac1 * (((1+z)/(1+this%z1))**(A01-A11+A21))*exp((A11-A21)*(z-this%z1)+A21*(z**2-this%z1**2)/2)
            else
                grho_de = grho_de_today * fac1 * fac2
            end if 
        else if(this%DEmodel == 10) then 
            ! Quadratic w(z): 3 bins  

            Delta_z1 = this%z1-0
            Delta_z2 = this%z2-this%z1
            Delta_z3 = this%z3-this%z2
            Delta_w1 = this%w1-this%w0
            Delta_w2 = this%w2-this%w1
            Delta_w3 = -1.0_dl-this%w2 

            wa0  = 2._dl*(Delta_w1/Delta_z1 - Delta_w2/Delta_z2 + Delta_w3/Delta_z3)
            waa0 = -Delta_w1/Delta_z1**2 + (2._dl/Delta_z1)*(Delta_w2/Delta_z2 - Delta_w3/Delta_z3)
            wa1  = 2._dl*(Delta_w2/Delta_z2 - Delta_w3/Delta_z3)
            waa1 = -Delta_w2/Delta_z2**2 + 2._dl*Delta_w3/(Delta_z2*Delta_z3)
            wa2  = 2._dl*Delta_w3/Delta_z3
            waa2 = -Delta_w3/Delta_z3**2

            A00 = 3._dl*(1 + this%w0 - wa0*0 + waa0*0**2._dl)
            A10 = 3._dl*(wa0 - 2*waa0*0)
            A20 = 3._dl*waa0

            A01 = 3._dl*(1 + this%w1 - wa1*this%z1 + waa1*this%z1**2._dl)
            A11 = 3._dl*(wa1 - 2*waa1*this%z1)
            A21 = 3._dl*waa1

            A02 = 3._dl*(1 + this%w2 - wa2*this%z2 + waa2*this%z2**2._dl)
            A12 = 3._dl*(wa2 - 2*waa2*this%z2)
            A22 = 3._dl*waa2

            fac1 = (((1+this%z1)/(1+0))**(A00-A10+A20))*exp((A10-A20)*(this%z1-0)+A20*(this%z1**2-0**2)/2)
            fac2 = (((1+this%z2)/(1+this%z1))**(A01-A11+A21))*exp((A11-A21)*(this%z2-this%z1)+A21*(this%z2**2-this%z1**2)/2)
            fac3 = (((1+this%z3)/(1+this%z2))**(A02-A12+A22))*exp((A12-A22)*(this%z3-this%z2)+A22*(this%z3**2-this%z2**2)/2)  

            if (z < this%z1) then
                grho_de = grho_de_today * (((1+z)/(1+0))**(A00-A10+A20))*exp((A10-A20)*(z -0)+A20*(z**2-0**2)/2)
            else if (z < this%z2) then
                grho_de = grho_de_today * fac1 * (((1+z)/(1+this%z1))**(A01-A11+A21))*exp((A11-A21)*(z-this%z1)+A21*(z**2-this%z1**2)/2)
            else if (z < this%z3) then
                grho_de = grho_de_today * fac1 * fac2 * (((1+z)/(1+this%z2))**(A02-A12+A22))*exp((A12-A22)*(z-this%z2)+A22*(z**2-this%z2**2)/2)
            else
                grho_de = grho_de_today * fac1 * fac2 * fac3
            end if     
        else if(this%DEmodel == 11) then 
            ! Quadratic w(z): 5 bins  

            Delta_z1 = this%z1-0
            Delta_z2 = this%z2-this%z1
            Delta_z3 = this%z3-this%z2
            Delta_z4 = this%z4-this%z3
            Delta_z5 = this%z5-this%z4

            Delta_w1 = this%w1-this%w0
            Delta_w2 = this%w2-this%w1
            Delta_w3 = this%w3-this%w2
            Delta_w4 = this%w4-this%w3
            Delta_w5 = -1.0_dl-this%w4 

            !See Mathematica Notebook - CAMBLateDE_boundary_conditions.nb
            wa0  = 2._dl*(Delta_w1/Delta_z1 - Delta_w2/Delta_z2 + Delta_w3/Delta_z3 - Delta_w4/Delta_z4 + Delta_w5/Delta_z5)
            waa0 = -Delta_w1/Delta_z1**2 + (2._dl/Delta_z1)*(Delta_w2/Delta_z2 - Delta_w3/Delta_z3 + Delta_w4/Delta_z4 - Delta_w5/Delta_z5)
            
            wa1  = 2._dl*(Delta_w2/Delta_z2 - Delta_w3/Delta_z3 + Delta_w4/Delta_z4 - Delta_w5/Delta_z5)
            waa1 = -Delta_w2/Delta_z2**2 + (2._dl/Delta_z2)*(Delta_w3/Delta_z3 - Delta_w4/Delta_z4 + Delta_w5/Delta_z5)
            
            wa2  = 2._dl*(Delta_w3/Delta_z3 - Delta_w4/Delta_z4 + Delta_w5/Delta_z5)
            waa2 = -Delta_w3/Delta_z3**2 + (2._dl/Delta_z3)*(Delta_w4/Delta_z4 - Delta_w5/Delta_z5)

            wa3  = 2._dl*(Delta_w4/Delta_z4 - Delta_w5/Delta_z5)
            waa3 = -Delta_w4/Delta_z4**2 + (2._dl/Delta_z4)*(Delta_w5/Delta_z5)

            wa4  = 2._dl*(Delta_w5/Delta_z5)
            waa4 = -Delta_w5/Delta_z5**2


            A00 = 3._dl*(1 + this%w0 - wa0*0 + waa0*0**2._dl)
            A10 = 3._dl*(wa0 - 2*waa0*0)
            A20 = 3._dl*waa0

            A01 = 3._dl*(1 + this%w1 - wa1*this%z1 + waa1*this%z1**2._dl)
            A11 = 3._dl*(wa1 - 2*waa1*this%z1)
            A21 = 3._dl*waa1

            A02 = 3._dl*(1 + this%w2 - wa2*this%z2 + waa2*this%z2**2._dl)
            A12 = 3._dl*(wa2 - 2*waa2*this%z2)
            A22 = 3._dl*waa2

            A03 = 3._dl*(1 + this%w3 - wa3*this%z3 + waa3*this%z3**2._dl)
            A13 = 3._dl*(wa3 - 2*waa3*this%z3)
            A23 = 3._dl*waa3

            A04 = 3._dl*(1 + this%w4 - wa4*this%z4 + waa4*this%z4**2._dl)
            A14 = 3._dl*(wa4 - 2*waa4*this%z4)
            A24 = 3._dl*waa4            

            fac1 = (((1+this%z1)/(1+0))**(A00-A10+A20))*exp((A10-A20)*(this%z1-0)+A20*(this%z1**2-0**2)/2)
            fac2 = (((1+this%z2)/(1+this%z1))**(A01-A11+A21))*exp((A11-A21)*(this%z2-this%z1)+A21*(this%z2**2-this%z1**2)/2)
            fac3 = (((1+this%z3)/(1+this%z2))**(A02-A12+A22))*exp((A12-A22)*(this%z3-this%z2)+A22*(this%z3**2-this%z2**2)/2)
            fac4 = (((1+this%z4)/(1+this%z3))**(A03-A13+A23))*exp((A13-A23)*(this%z4-this%z3)+A23*(this%z4**2-this%z3**2)/2)
            fac5 = (((1+this%z5)/(1+this%z4))**(A04-A14+A24))*exp((A14-A24)*(this%z5-this%z4)+A24*(this%z5**2-this%z4**2)/2)  

            if (z < this%z1) then
                grho_de = grho_de_today * (((1+z)/(1+0))**(A00-A10+A20))*exp((A10-A20)*(z -0)+A20*(z**2-0**2)/2)
            else if (z < this%z2) then
                grho_de = grho_de_today * fac1 * (((1+z)/(1+this%z1))**(A01-A11+A21))*exp((A11-A21)*(z-this%z1)+A21*(z**2-this%z1**2)/2)
            else if (z < this%z3) then
                grho_de = grho_de_today * fac1 * fac2 * (((1+z)/(1+this%z2))**(A02-A12+A22))*exp((A12-A22)*(z-this%z2)+A22*(z**2-this%z2**2)/2)
            else if (z < this%z4) then
                grho_de = grho_de_today * fac1 * fac2 * fac3 * (((1+z)/(1+this%z3))**(A03-A13+A23))*exp((A13-A23)*(z-this%z3)+A23*(z**2-this%z3**2)/2)
            else if (z < this%z5) then    
                grho_de = grho_de_today * fac1 * fac2 * fac3 * fac4 * (((1+z)/(1+this%z4))**(A04-A14+A24))*exp((A14-A24)*(z-this%z4)+A24*(z**2-this%z4**2)/2)
            else
                grho_de = grho_de_today * fac1 * fac2 * fac3 * fac4 * fac5
            end if                     
        else if(this%DEmodel == 12) then 
            ! Cubic w(z): 2 bins

            Delta_z1 = this%z1-0
            Delta_z2 = this%z2-this%z1
            Delta_w1 = this%w1-this%w0
            Delta_w2 = -1.0_dl-this%w1

            wa0   = 3*Delta_w1/Delta_z1 - 3*Delta_w2/Delta_z2 * (2+Delta_z1/Delta_z2)
            waa0  = -3*Delta_w1/Delta_z1**2 + 3*Delta_w2/Delta_z2*(3/Delta_z1 + 2/Delta_z2)
            waaa0 = Delta_w1/Delta_z1**3 - 3*Delta_w2/(Delta_z1*Delta_z2) * (1/Delta_z1 + 1/Delta_z2)
            wa1   = 3*Delta_w2/Delta_z2
            waa1  = -3*Delta_w2/Delta_z2**2
            waaa1 = Delta_w2/Delta_z2**3

            A00 = 3*(1 + this%w0 - wa0*0 + waa0*0**2 - waaa0*0**3)
            A10 = 3*(wa0 - 2*waa0*0 + 3*waaa0*0**2)
            A20 = 3*(waa0 - 3*waaa0*0)
            A30 = 3*waaa0

            A01 = 3*(1 + this%w1 - wa1*this%z1 + waa1*this%z1**2 - waaa1*this%z1**3)
            A11 = 3*(wa1 - 2*waa1*this%z1 + 3*waaa1*this%z1**2)
            A21 = 3*(waa1 - 3*waaa1*this%z1)
            A31 = 3*waaa1

            fac1 = (((1+this%z1)/(1+0))**(A00-A10+A20-A30)) * &
                   exp((A10-A20+A30)*(this%z1-0)+(A20-A30)*(this%z1**2-0**2)/2+A30*(this%z1**3-0**3)/3)

            fac2 = (((1+this%z2)/(1+this%z1))**(A01-A11+A21-A31)) * &
                   exp((A11-A21+A31)*(this%z2-this%z1)+(A21-A31)*(this%z2**2-this%z1**2)/2+A31*(this%z2**3-this%z1**3)/3)

            if (z < this%z1) then
                grho_de = grho_de_today * (((1+z)/(1+0))**(A00-A10+A20-A30)) * &
                          exp((A10-A20+A30)*(z-0)+(A20-A30)*(z**2-0**2)/2+A30*(z**3-0**3)/3)
            else if (z < this%z2) then
                grho_de = grho_de_today * fac1 * (((1+z)/(1+this%z1))**(A01-A11+A21-A31)) * &
                          exp((A11-A21+A31)*(z-this%z1)+(A21-A31)*(z**2-this%z1**2)/2+A31*(z**3-this%z1**3)/3)
            else
                grho_de = grho_de_today * fac1 * fac2
            end if
        else if(this%DEmodel == 13) then 
            ! Cubic w(z): 3 bins

            Delta_z1 = this%z1-0
            Delta_z2 = this%z2-this%z1
            Delta_z3 = this%z3-this%z2
            Delta_w1 = this%w1-this%w0
            Delta_w2 = this%w2-this%w1
            Delta_w3 = -1.0_dl-this%w2 

            wa0   = 3*Delta_w1/Delta_z1 - 3*(Delta_w2/Delta_z2)*(2 + Delta_z1/Delta_z2) + 3*(Delta_w3/Delta_z3)*(4 + 3*Delta_z1/Delta_z2 + 2*Delta_z1/Delta_z3 + 2*Delta_z2/Delta_z3)
            waa0  = -3*Delta_w1/Delta_z1**2 + 3*(Delta_w2/Delta_z2)*(3/Delta_z1 + 2/Delta_z2) -3*(Delta_w3/Delta_z3)*(6/Delta_z1 + 6/Delta_z2 + 4/Delta_z3 + 3*Delta_z2/(Delta_z1*Delta_z3))
            waaa0 = Delta_w1/Delta_z3**3 - 3*Delta_w2/(Delta_z1*Delta_z2)*(1/Delta_z1 + 1/Delta_z2) + 3*Delta_w3/(Delta_z1*Delta_z3)*(2/Delta_z1 + 3/Delta_z2 + 2/Delta_z3 + Delta_z2/(Delta_z1*Delta_z3))
            wa1   = 3*Delta_w2/Delta_z2 - 3*Delta_w3/Delta_z3*(2 + Delta_z2/Delta_z3)
            waa1  = -3*Delta_w2/Delta_z2**2 + 3*Delta_w3/Delta_z3*(3/Delta_z2 + 2/Delta_z3)
            waaa1 = Delta_w2/Delta_z2**3 - 3*Delta_w3/(Delta_z2*Delta_z3)*(1/Delta_z2 + 1/Delta_z3)
            wa2   = 3*Delta_w3/Delta_z3
            waa2  = -3*Delta_w3/Delta_z3**2
            waaa2 = Delta_w3/Delta_z3**3

            A00 = 3*(1 + this%w0 - wa0*0 + waa0*0**2 - waaa0*0**3)
            A10 = 3*(wa0 - 2*waa0*0 + 3*waaa0*0**2)
            A20 = 3*(waa0 - 3*waaa0*0)
            A30 = 3*waaa0

            A01 = 3*(1 + this%w1 - wa1*this%z1 + waa1*this%z1**2 - waaa1*this%z1**3)
            A11 = 3*(wa1 - 2*waa1*this%z1 + 3*waaa1*this%z1**2)
            A21 = 3*(waa1 - 3*waaa1*this%z1)
            A31 = 3*waaa1

            A02 = 3*(1 + this%w2 - wa2*this%z2 + waa2*this%z2**2 - waaa2*this%z2**3)
            A12 = 3*(wa2 - 2*waa2*this%z2 + 3*waaa2*this%z2**2)
            A22 = 3*(waa2 - 3*waaa2*this%z2)
            A32 = 3*waaa2

            fac1 = (((1+this%z1)/(1+0))**(A00-A10+A20-A30)) * &
                   exp((A10-A20+A30)*(this%z1-0)+(A20-A30)*(this%z1**2-0**2)/2+A30*(this%z1**3-0**3)/3)

            fac2 = (((1+this%z2)/(1+this%z1))**(A01-A11+A21-A31)) * &
                   exp((A11-A21+A31)*(this%z2-this%z1)+(A21-A31)*(this%z2**2-this%z1**2)/2+A31*(this%z2**3-this%z1**3)/3)

            fac3 = (((1+this%z3)/(1+this%z2))**(A02-A12+A22-A32)) * &
                   exp((A12-A22+A32)*(this%z3-this%z2)+(A22-A32)*(this%z3**2-this%z2**2)/2+A32*(this%z3**3-this%z2**3)/3)

            if (z < this%z1) then
                grho_de = grho_de_today * (((1+z)/(1+0))**(A00-A10+A20-A30)) * &
                          exp((A10-A20+A30)*(z-0)+(A20-A30)*(z**2-0**2)/2+A30*(z**3-0**3)/3)
            else if (z < this%z2) then
                grho_de = grho_de_today * fac1 * (((1+z)/(1+this%z1))**(A01-A11+A21-A31)) * &
                          exp((A11-A21+A31)*(z-this%z1)+(A21-A31)*(z**2-this%z1**2)/2+A31*(z**3-this%z1**3)/3)
            else if (z < this%z3) then
                grho_de = grho_de_today * fac1 * fac2 * (((1+z)/(1+this%z2))**(A02-A12+A22-A32)) * &
                          exp((A12-A22+A32)*(z-this%z2)+(A22-A32)*(z**2-this%z2**2)/2+A32*(z**3-this%z2**3)/3)
            else
                grho_de = grho_de_today * fac1 * fac2 * fac3                           
            end if

        else if(this%DEmodel == 14) then 
            ! Cubic w(z): 5 bins

            Delta_z1 = this%z1-0
            Delta_z2 = this%z2-this%z1
            Delta_z3 = this%z3-this%z2
            Delta_z4 = this%z4-this%z3
            Delta_z5 = this%z5-this%z4

            Delta_w1 = this%w1-this%w0
            Delta_w2 = this%w2-this%w1
            Delta_w3 = this%w3-this%w2
            Delta_w4 = this%w4-this%w3
            Delta_w5 = -1.0_dl-this%w4 

            !See Mathematica Notebook - CAMBLateDE_boundary_conditions.nb
            wa0   = (3*(Delta_w1))/(Delta_z1) - (3*(Delta_w2)*(Delta_z1))/(Delta_z2)**2 - (6*(Delta_w2))/(Delta_z2) + &
                    ((Delta_w3)*(6*(Delta_z1) + 6*(Delta_z2)))/(Delta_z3)**2 + ((Delta_w3)*(12 + (9*(Delta_z1))/(Delta_z2)))/(Delta_z3) + &
                    ((Delta_w4)*(-12*(Delta_z1) - 12*(Delta_z2) - 12*(Delta_z3) - (9*(Delta_z1)*(Delta_z3))/(Delta_z2)))/(Delta_z4)**2 + &
                    ((Delta_w4)*(-24 - (18*(Delta_z1))/(Delta_z2) - (18*(Delta_z1))/(Delta_z3) - (18*(Delta_z2))/(Delta_z3)))/(Delta_z4) + & 
                    ((Delta_w5)*(24*(Delta_z1) + 24*(Delta_z2) + 24*(Delta_z3) + (18*(Delta_z1)*(Delta_z3))/(Delta_z2) + 24*(Delta_z4) + &
                    (18*(Delta_z1)*(Delta_z4))/(Delta_z2) + (18*(Delta_z1)*(Delta_z4))/(Delta_z3) + (18*(Delta_z2)*(Delta_z4))/(Delta_z3)))/(Delta_z5)**2 + &
                    ((Delta_w5)*(48 + (36*(Delta_z1))/(Delta_z2) + (36*(Delta_z1))/(Delta_z3) + (36*(Delta_z2))/(Delta_z3) + (36*(Delta_z1))/(Delta_z4) + &
                    (36*(Delta_z2))/(Delta_z4) + (36*(Delta_z3))/(Delta_z4) + (27*(Delta_z1)*(Delta_z3))/((Delta_z2)*(Delta_z4))))/(Delta_z5)
            waa0  = (-3*(Delta_w1))/(Delta_z1)**2 + (6*(Delta_w2))/(Delta_z2)**2 + (9*(Delta_w2))/((Delta_z1)*(Delta_z2)) + &
                    ((Delta_w3)*(-12 - (9*(Delta_z2))/(Delta_z1)))/(Delta_z3)**2 + ((Delta_w3)*(-18/(Delta_z1) - 18/(Delta_z2)))/(Delta_z3) + &
                    ((Delta_w4)*(24 + (18*(Delta_z2))/(Delta_z1) + (18*(Delta_z3))/(Delta_z1) + (18*(Delta_z3))/(Delta_z2)))/(Delta_z4)**2 + &
                    ((Delta_w4)*(36/(Delta_z1) + 36/(Delta_z2) + 36/(Delta_z3) + (27*(Delta_z2))/((Delta_z1)*(Delta_z3))))/(Delta_z4) + &
                    ((Delta_w5)*(-48 - (36*(Delta_z2))/(Delta_z1) - (36*(Delta_z3))/(Delta_z1) - (36*(Delta_z3))/(Delta_z2) - (36*(Delta_z4))/(Delta_z1) - &
                    (36*(Delta_z4))/(Delta_z2) - (36*(Delta_z4))/(Delta_z3) - (27*(Delta_z2)*(Delta_z4))/((Delta_z1)*(Delta_z3))))/(Delta_z5)**2 + &
                    ((Delta_w5)*(-72/(Delta_z1) - 72/(Delta_z2) - 72/(Delta_z3) - (54*(Delta_z2))/((Delta_z1)*(Delta_z3)) - 72/(Delta_z4) - &
                    (54*(Delta_z2))/((Delta_z1)*(Delta_z4)) - (54*(Delta_z3))/((Delta_z1)*(Delta_z4)) - (54*(Delta_z3))/((Delta_z2)*(Delta_z4))))/(Delta_z5)
            waaa0 = (Delta_w1)/(Delta_z1)**3 - (3*(Delta_w2))/((Delta_z1)*(Delta_z2)**2) - (3*(Delta_w2))/((Delta_z1)**2*(Delta_z2)) + &
                    ((Delta_w3)*(6/(Delta_z1) + (3*(Delta_z2))/(Delta_z1)**2))/(Delta_z3)**2 + ((Delta_w3)*(6/(Delta_z1)**2 + 9/((Delta_z1)*(Delta_z2))))/(Delta_z3) + &
                    ((Delta_w4)*(-12/(Delta_z1) - (6*(Delta_z2))/(Delta_z1)**2 - (6*(Delta_z3))/(Delta_z1)**2 - (9*(Delta_z3))/((Delta_z1)*(Delta_z2))))/&
                    (Delta_z4)**2 + ((Delta_w4)*(-12/(Delta_z1)**2 - 18/((Delta_z1)*(Delta_z2)) - 18/((Delta_z1)*(Delta_z3)) - &
                    (9*(Delta_z2))/((Delta_z1)**2*(Delta_z3))))/(Delta_z4) + &
                    ((Delta_w5)*(24/(Delta_z1) + (12*(Delta_z2))/(Delta_z1)**2 + (12*(Delta_z3))/(Delta_z1)**2 + (18*(Delta_z3))/((Delta_z1)*(Delta_z2)) + &
                    (12*(Delta_z4))/(Delta_z1)**2 + (18*(Delta_z4))/((Delta_z1)*(Delta_z2)) + (18*(Delta_z4))/((Delta_z1)*(Delta_z3)) + &
                    (9*(Delta_z2)*(Delta_z4))/((Delta_z1)**2*(Delta_z3))))/(Delta_z5)**2 + &
                    ((Delta_w5)*(24/(Delta_z1)**2 + 36/((Delta_z1)*(Delta_z2)) + 36/((Delta_z1)*(Delta_z3)) + (18*(Delta_z2))/((Delta_z1)**2*(Delta_z3)) + &
                    36/((Delta_z1)*(Delta_z4)) + (18*(Delta_z2))/((Delta_z1)**2*(Delta_z4)) + (18*(Delta_z3))/((Delta_z1)**2*(Delta_z4)) + &
                    (27*(Delta_z3))/((Delta_z1)*(Delta_z2)*(Delta_z4))))/(Delta_z5)
            
            wa1   = (3*(Delta_w2))/(Delta_z2) - (3*(Delta_w3)*(Delta_z2))/(Delta_z3)**2 - (6*(Delta_w3))/(Delta_z3) + &
                    ((Delta_w4)*(6*(Delta_z2) + 6*(Delta_z3)))/(Delta_z4)**2 + ((Delta_w4)*(12 + (9*(Delta_z2))/(Delta_z3)))/(Delta_z4) + &
                    ((Delta_w5)*(-12*(Delta_z2) - 12*(Delta_z3) - 12*(Delta_z4) - (9*(Delta_z2)*(Delta_z4))/(Delta_z3)))/(Delta_z5)**2 + &
                    ((Delta_w5)*(-24 - (18*(Delta_z2))/(Delta_z3) - (18*(Delta_z2))/(Delta_z4) - (18*(Delta_z3))/(Delta_z4)))/(Delta_z5)
            waa1  = (-3*(Delta_w2))/(Delta_z2)**2 + (6*(Delta_w3))/(Delta_z3)**2 + (9*(Delta_w3))/((Delta_z2)*(Delta_z3)) + &
                    ((Delta_w4)*(-12 - (9*(Delta_z3))/(Delta_z2)))/(Delta_z4)**2 + ((Delta_w4)*(-18/(Delta_z2) - 18/(Delta_z3)))/(Delta_z4) + &
                    ((Delta_w5)*(24 + (18*(Delta_z3))/(Delta_z2) + (18*(Delta_z4))/(Delta_z2) + (18*(Delta_z4))/(Delta_z3)))/(Delta_z5)**2 + &
                    ((Delta_w5)*(36/(Delta_z2) + 36/(Delta_z3) + 36/(Delta_z4) + (27*(Delta_z3))/((Delta_z2)*(Delta_z4))))/(Delta_z5)
            waaa1 = (Delta_w2)/(Delta_z2)**3 - (3*(Delta_w3))/((Delta_z2)*(Delta_z3)**2) - (3*(Delta_w3))/((Delta_z2)**2*(Delta_z3)) + &
                    ((Delta_w4)*(6/(Delta_z2) + (3*(Delta_z3))/(Delta_z2)**2))/(Delta_z4)**2 + ((Delta_w4)*(6/(Delta_z2)**2 + 9/((Delta_z2)*(Delta_z3))))/(Delta_z4) + &
                    ((Delta_w5)*(-12/(Delta_z2) - (6*(Delta_z3))/(Delta_z2)**2 - (6*(Delta_z4))/(Delta_z2)**2 - (9*(Delta_z4))/((Delta_z2)*(Delta_z3))))/&
                    (Delta_z5)**2 + ((Delta_w5)*(-12/(Delta_z2)**2 - 18/((Delta_z2)*(Delta_z3)) - 18/((Delta_z2)*(Delta_z4)) - &
                    (9*(Delta_z3))/((Delta_z2)**2*(Delta_z4))))/(Delta_z5)
            
            wa2   = (3*(Delta_w3))/(Delta_z3) - (3*(Delta_w4)*(Delta_z3))/(Delta_z4)**2 - (6*(Delta_w4))/(Delta_z4) + &
                    ((Delta_w5)*(6*(Delta_z3) + 6*(Delta_z4)))/(Delta_z5)**2 + ((Delta_w5)*(12 + (9*(Delta_z3))/(Delta_z4)))/(Delta_z5)
            waa2  = (-3*(Delta_w3))/(Delta_z3)**2 + (6*(Delta_w4))/(Delta_z4)**2 + (9*(Delta_w4))/((Delta_z3)*(Delta_z4)) + &
                    ((Delta_w5)*(-12 - (9*(Delta_z4))/(Delta_z3)))/(Delta_z5)**2 + ((Delta_w5)*(-18/(Delta_z3) - 18/(Delta_z4)))/(Delta_z5)
            waaa2 = (Delta_w3)/(Delta_z3)**3 - (3*(Delta_w4))/((Delta_z3)*(Delta_z4)**2) - (3*(Delta_w4))/((Delta_z3)**2*(Delta_z4)) + &
                    ((Delta_w5)*(6/(Delta_z3) + (3*(Delta_z4))/(Delta_z3)**2))/(Delta_z5)**2 + ((Delta_w5)*(6/(Delta_z3)**2 + 9/((Delta_z3)*(Delta_z4))))/(Delta_z5)

            wa3   = (3*(Delta_w4))/(Delta_z4) - (3*(Delta_w5)*(Delta_z4))/(Delta_z5)**2 - (6*(Delta_w5))/(Delta_z5)
            waa3  = (-3*(Delta_w4))/(Delta_z4)**2 + (6*(Delta_w5))/(Delta_z5)**2 + (9*(Delta_w5))/((Delta_z4)*(Delta_z5))
            waaa3 = (Delta_w4)/(Delta_z4)**3 - (3*(Delta_w5))/((Delta_z4)*(Delta_z5)**2) - (3*(Delta_w5))/((Delta_z4)**2*(Delta_z5))

            wa4   = (3*(Delta_w5))/(Delta_z5)
            waa4  = (-3*(Delta_w5))/(Delta_z5)**2
            waaa4 = (Delta_w5)/(Delta_z5)**3

            A00 = 3*(1 + this%w0 - wa0*0 + waa0*0**2 - waaa0*0**3)
            A10 = 3*(wa0 - 2*waa0*0 + 3*waaa0*0**2)
            A20 = 3*(waa0 - 3*waaa0*0)
            A30 = 3*waaa0

            A01 = 3*(1 + this%w1 - wa1*this%z1 + waa1*this%z1**2 - waaa1*this%z1**3)
            A11 = 3*(wa1 - 2*waa1*this%z1 + 3*waaa1*this%z1**2)
            A21 = 3*(waa1 - 3*waaa1*this%z1)
            A31 = 3*waaa1

            A02 = 3*(1 + this%w2 - wa2*this%z2 + waa2*this%z2**2 - waaa2*this%z2**3)
            A12 = 3*(wa2 - 2*waa2*this%z2 + 3*waaa2*this%z2**2)
            A22 = 3*(waa2 - 3*waaa2*this%z2)
            A32 = 3*waaa2

            A03 = 3*(1 + this%w3 - wa3*this%z3 + waa3*this%z3**2 - waaa3*this%z3**3)
            A13 = 3*(wa3 - 2*waa3*this%z3 + 3*waaa3*this%z3**2)
            A23 = 3*(waa3 - 3*waaa3*this%z3)
            A33 = 3*waaa3

            A04 = 3*(1 + this%w4 - wa4*this%z4 + waa4*this%z4**2 - waaa4*this%z4**3)
            A14 = 3*(wa4 - 2*waa4*this%z4 + 3*waaa4*this%z4**2)
            A24 = 3*(waa4 - 3*waaa4*this%z4)
            A34 = 3*waaa4             

            fac1 = (((1+this%z1)/(1+0))**(A00-A10+A20-A30)) * &
                   exp((A10-A20+A30)*(this%z1-0)+(A20-A30)*(this%z1**2-0**2)/2+A30*(this%z1**3-0**3)/3)

            fac2 = (((1+this%z2)/(1+this%z1))**(A01-A11+A21-A31)) * &
                   exp((A11-A21+A31)*(this%z2-this%z1)+(A21-A31)*(this%z2**2-this%z1**2)/2+A31*(this%z2**3-this%z1**3)/3)

            fac3 = (((1+this%z3)/(1+this%z2))**(A02-A12+A22-A32)) * &
                   exp((A12-A22+A32)*(this%z3-this%z2)+(A22-A32)*(this%z3**2-this%z2**2)/2+A32*(this%z3**3-this%z2**3)/3)

            fac4 = (((1+this%z4)/(1+this%z3))**(A03-A13+A23-A33)) * &
                   exp((A13-A23+A33)*(this%z4-this%z3)+(A23-A33)*(this%z4**2-this%z3**2)/2+A33*(this%z4**3-this%z3**3)/3)

            fac5 = (((1+this%z5)/(1+this%z4))**(A04-A14+A24-A34)) * &
                   exp((A14-A24+A34)*(this%z5-this%z4)+(A24-A34)*(this%z5**2-this%z4**2)/2+A34*(this%z5**3-this%z4**3)/3)

            if (z < this%z1) then
                grho_de = grho_de_today * & 
                          (((1+z)/(1+0))**(A00-A10+A20-A30)) * &
                          exp((A10-A20+A30)*(z-0)+(A20-A30)*(z**2-0**2)/2+A30*(z**3-0**3)/3)
            else if (z < this%z2) then
                grho_de = grho_de_today * fac1 * & 
                          (((1+z)/(1+this%z1))**(A01-A11+A21-A31)) * &
                          exp((A11-A21+A31)*(z-this%z1)+(A21-A31)*(z**2-this%z1**2)/2+A31*(z**3-this%z1**3)/3)
            else if (z < this%z3) then
                grho_de = grho_de_today * fac1 * fac2 * & 
                          (((1+z)/(1+this%z2))**(A02-A12+A22-A32)) * &
                          exp((A12-A22+A32)*(z-this%z2)+(A22-A32)*(z**2-this%z2**2)/2+A32*(z**3-this%z2**3)/3)
            else if (z < this%z4) then
                grho_de = grho_de_today * fac1 * fac2 * fac3 * &
                          (((1+z)/(1+this%z3))**(A03-A13+A23-A33)) * &
                          exp((A13-A23+A33)*(z-this%z3)+(A23-A33)*(z**2-this%z3**2)/2+A33*(z**3-this%z3**3)/3)            
            else if (z < this%z5) then
                grho_de = grho_de_today * fac1 * fac2 * fac3 * fac4 * & 
                          (((1+z)/(1+this%z4))**(A04-A14+A24-A34)) * &
                          exp((A14-A24+A34)*(z-this%z4)+(A24-A34)*(z**2-this%z4**2)/2+A34*(z**3-this%z4**3)/3)             
            else
                grho_de = grho_de_today * fac1 * fac2 * fac3 * fac4 * fac5                           
            end if            
        else if (this%DEmodel == 15) then 
            ! tanh model
            a1 = (this%w1 - this%w0)/this%sigma/2._dl
            a2 = this%w1 - a1*(this%z1 + this%sigma)
            za = this%z1 - this%sigma
            zb = this%z1 + this%sigma
            if (z < this%z1 - this%sigma) then
                grho_de = grho_de_today * (1.0_dl+z)**(3._dl*(1._dl + this%w0))
            else if (z > this%z1 + this%sigma) then
                grho_de = grho_de_today * (1.0_dl + za)**(3._dl*(1+this%w0)) * &
                 exp(a1*(zb-za)) * ((1.0_dl + zb)/(1.0_dl + za))**(1.0_dl - a1 + a2) * &
                 ((1.0_dl+z)/(1.0_dl + zb))**(3.0_dl*(1.0_dl + this%w1))
            else
                grho_de = grho_de_today * (1.0_dl + za)**(3._dl*(1+this%w0)) * &
                 ( exp(a1*(z-za)) * ((1.0_dl + z)/(1.0_dl + za))**(1.0_dl-a1+a2))                
            end if
        !DHFS MOD TANH START    
        else if (this%DEmodel == 16) then   
            ! Numeric tanh         
            if (z < this%z1+5.0*this%sigma) then
                grho_de = grho_de_today * exp( 3.0_dl * Integrate_Romberg(this, kernel_tanh, 0.0_dl, z, 1d-5, 20, 100) ) 
            else 
                grho_de = grho_de_today *(z/(this%z1+5.0_dl*this%sigma))**(3.0_dl*(1.0_dl+this%w1)) * exp( 3.0_dl * Integrate_tanh ) 
            end if                    
        !DHFS MOD TANH START    
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
        !DHFS MOD TANH START
        if (this%DEmodel == 16) then
            ! Numeric tanh
            Integrate_tanh = Integrate_Romberg(this, kernel_tanh, 0.0_dl, this%z1+5.0*this%sigma, 1d-5, 20, 100)
        end if
        !DHFS MOD TANH END
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

        if (this%DEmodel == 2) then
            !'w0wa'
            w  = this%w0
            wa = this%w1
        else if (this%DEmodel==1  .or. this%DEmodel==3  .or. this%DEmodel==4  .or. this%DEmodel==5  .or. &
                 this%DEmodel==6  .or. this%DEmodel==7  .or. this%DEmodel==8  .or. this%DEmodel==9  .or. &
                 this%DEmodel==10 .or. this%DEmodel==11 .or. this%DEmodel==12 .or. this%DEmodel==13 .or. & 
                 this%DEmodel==14) then
            w  = this%w0
            wa = 0
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
