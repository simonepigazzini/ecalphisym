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

    def sum(self, axis: int=0):
        """
        Sum an array of Info

        :param axis: axis along which the sum is performed.
        """
        data = {
            "hitseb": ak.sum(self.hitseb, axis=axis),
            "hitsee": ak.sum(self.hitsee, axis=axis),
            "nevents": ak.sum(self.nevents, axis=axis),
            "nlumis": ak.sum(self.nlumis, axis=axis),
            "fill": np.unique(self.fill),
            "reclumi": ak.sum(self.reclumi, axis=axis),
            "delivlumi": ak.sum(self.delivlumi, axis=axis)
        }
        return ak.zip(
            data,
            with_name="Info"
        )
    
__all__ = ["Info"]
