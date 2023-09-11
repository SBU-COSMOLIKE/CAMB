from .baseconfig import F2003Class, fortran_class, numpy_1d, CAMBError, np, \
    AllocatableArrayDouble, f_pointer
from ctypes import c_int, c_double, byref, POINTER, c_bool


class DarkEnergyModel(F2003Class):
    """
    Abstract base class for dark energy model implementations.
    """
    _fields_ = [
        ("__is_cosmological_constant", c_bool),
        ("__num_perturb_equations", c_int)]

    def validate_params(self):
        return True

@fortran_class
class LateDE(DarkEnergyModel):
    """
    DHFS: Abstract base class for models using w and wa parameterization with use w(a) = w + (1-a)*wa parameterization,
    or bin w.
    """

    _fortran_class_module_ = 'LateDE'
    _fortran_class_name_ = 'TLateDE'

    _fields_ = [
        ("model", c_int, "select one model among five: (1) w=cte, (2) CPL, (3) 3 bins w, (4) 5 bins w (5) 10 bins w"),
        # Equation of State
        ("w0", c_double, "Bin w parameter: EoS for the 0th bin"),
        ("w1", c_double, "Bin w parameter: EoS for the 1st bin"),
        ("w2", c_double, "Bin w parameter: EoS for the 2nd bin"),
        ("w3", c_double, "Bin w parameter: EoS for the 3rd bin"),
        ("w4", c_double, "Bin w parameter: EoS for the 4th bin"),
        ("w5", c_double, "Bin w parameter: EoS for the 5th bin"),
        ("w6", c_double, "Bin w parameter: EoS for the 6th bin"),
        ("w7", c_double, "Bin w parameter: EoS for the 7th bin"),
        ("w8", c_double, "Bin w parameter: EoS for the 8th bin"),
        ("w9", c_double, "Bin w parameter: EoS for the 9th bin"),
        ("w10", c_double, "Bin w parameter: EoS for the 10th bin"),
        # Redshift 
        ("z1", c_double, "Bin w parameter: redshift for the 1st bin"),
        ("z2", c_double, "Bin w parameter: redshift for the 2nd bin"),
        ("z3", c_double, "Bin w parameter: redshift for the 3rd bin"),
        ("z4", c_double, "Bin w parameter: redshift for the 4th bin"),
        ("z5", c_double, "Bin w parameter: redshift for the 5th bin"),
        ("z6", c_double, "Bin w parameter: redshift for the 6th bin"),
        ("z7", c_double, "Bin w parameter: redshift for the 7th bin"),
        ("z8", c_double, "Bin w parameter: redshift for the 8th bin"),
        ("z9", c_double, "Bin w parameter: redshift for the 9th bin"),
        ("z10", c_double, " Bin w parameter: redshift for the 10th bin"),
        # Factors to match the integrated boundary condition for energy density
        ("fac1", c_double, "Bin w internal parameter: integrated boundary condition between the 1st and 2nd bin"),
        ("fac2", c_double, "Bin w internal parameter: integrated boundary condition between the 2st and 2nd bin"),
        ("fac3", c_double, "Bin w internal parameter: integrated boundary condition between the 3st and 2nd bin"),
        ("fac4", c_double, "Bin w internal parameter: integrated boundary condition between the 4st and 2nd bin"),
        ("fac5", c_double, "Bin w internal parameter: integrated boundary condition between the 5st and 2nd bin"),
        ("fac6", c_double, "Bin w internal parameter: integrated boundary condition between the 6st and 2nd bin"),
        ("fac7", c_double, "Bin w internal parameter: integrated boundary condition between the 7st and 2nd bin"),
        ("fac8", c_double, "Bin w internal parameter: integrated boundary condition between the 8st and 2nd bin"),
        ("fac9", c_double, "Bin w internal parameter: integrated boundary condition between the 9st and 2nd bin"),
        ("fac10", c_double, "Bin w internal parameter: integrated boundary condition between the 1st and 2nd bin")
    ]

    def set_params(self, model,
                     w0=-1, w1=-1, w2=-1, w3=-1, w4=-1, w5=-1, w6=-1, w7=-1, w8=-1, w9=-1, w10=-1,
                     z1=1, z2=2, z3=3, z4=4, z5=5, z6=6, z7=7, z8=8, z9=9, z10=10):

        self.model=model 
        self.w0=w0 
        self.w1=w1 
        self.w2=w2
        self.w3=w3
        self.w4=w4
        self.w5=w5
        self.w6=w6
        self.w7=w7
        self.w8=w8
        self.w9=w9
        self.w10=w10
        self.z1=z1
        self.z2=z2 
        self.z3=z3
        self.z4=z4 
        self.z5=z5
        self.z6=z6 
        self.z7=z7
        self.z8=z8
        self.z9=z9
        self.z10=z10

@fortran_class
class DarkEnergyPPF(LateDE):
    """
    VM: CLASS IMPLEMENTS w, w0wa and binw

    """
    # cannot declare c_Gamma_ppf directly here as have not defined all fields in DarkEnergyEqnOfState (TCubicSpline)
    _fortran_class_module_ = 'DarkEnergyPPF'
    _fortran_class_name_ = 'TDarkEnergyPPF'


# short names for models that support w/wa
F2003Class._class_names.update({'ppf': DarkEnergyPPF})
