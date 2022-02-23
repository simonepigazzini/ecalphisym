"""EcalPhiSymInfo with numpy mixins

This provides a class to handle an EcalPhiSymInfo and
implements high level operations among its data members
"""
import numpy as np
import awkward as ak

behavior = {}

@ak.mixin_class(behavior)
class Info:
    """
    An EcalPhiSymInfo contains global information concerning a Run/LuminosityBlock.
    Both LHC and ECAL related global quantities are stored in EcalPhiSymInfo.
    """

    @property
    def miscalibs_eb(self):
        """
        An array of the applied mis-calibration values for EB.
        The nominal value corrispond to zero (len/2-th element)
        """
        # linspace requires an integer singleton
        # while nmis is an array (one int per Run/LuminosityBlock        
        return np.linspace(ak.broadcast_arrays(self.minmiseb, np.ones(61200))[0]-1,
                           ak.broadcast_arrays(self.maxmiseb, np.ones(61200))[0]-1,
                           int(ak.mean(self.nmis))+1,
                           axis=-1)
    
    @property
    def miscalibs_ee(self):
        """
        An array of the applied mis-calibration values for EE.
        The nominal value corrispond to zero (len/2-th element)
        """
        # linspace requires an integer singleton
        # while nmis is an array (one int per Run/LuminosityBlock
        return np.linspace(ak.flatten(self.minmisee),
                           ak.flatten(self.maxmisee),
                           int(ak.mean(self.nmis))+1,
                           axis=-1)
    
__all__ = ["Info"]
