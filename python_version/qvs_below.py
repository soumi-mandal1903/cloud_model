import numpy as np

class QVSBelowCalculator:
    def __init__(self, qv_dtdlnp, qv_p, qv_t, qv_factor, qv_gas_name, vapor_pressure_obj):
        """
        vapor_pressure_obj: an object with a pvap_gas(gas_name, t, p) method
        """
        self.qv_dtdlnp = qv_dtdlnp
        self.qv_p = qv_p
        self.qv_t = qv_t
        self.qv_factor = qv_factor
        self.qv_gas_name = qv_gas_name
        self.vapor_pressure_obj = vapor_pressure_obj

    def qvs_below(self, p_test):
        """
        Calculate saturation mixing ratio for a gas, extrapolated below model domain.
        Input: p_test (pressure in dyne/cm^2)
        Output: saturation mixing ratio
        """
        t_test = self.qv_t + np.log(self.qv_p / p_test) * self.qv_dtdlnp
        pvap_test = self.vapor_pressure_obj.pvap_gas(self.qv_gas_name, t_test, p_test)
        qvs_below = self.qv_factor * pvap_test / p_test
        return qvs_below

# Example usage:
# vapor_pressure = VaporPressure()  # Your class with pvap_gas implemented
# qvs_calc = QVSBelowCalculator(qv_dtdlnp, qv_p, qv_t, qv_factor, qv_gas_name, vapor_pressure)
# result = qvs_calc.qvs_below(p_test)