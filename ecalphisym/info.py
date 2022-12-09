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
        return np.linspace(ak.broadcast_arrays(self.minmiseb.to_numpy(), np.ones(61200))[0]-1,
                           ak.broadcast_arrays(self.maxmiseb.to_numpy(), np.ones(61200))[0]-1,
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
        Sum an array of Info. The summing operation assumes that the miscalibration
        values are the same across the analyzed dataset. Further more if summing across
        fills the sum fill field is filled with the number of the lowest fill.

        :param axis: axis along which the sum is performed.
        """
        data = {
            'maxmiseb' : ak.mean(self.maxmiseb, axis=axis),
            'maxmisee' : ak.mean(self.maxmisee, axis=axis),
            'minmiseb': ak.mean(self.minmiseb, axis=axis),
            'minmisee': ak.mean(self.minmisee, axis=axis),
            'nmis' : ak.mean(self.nmis, axis=axis),
            'hitseb': ak.sum(self.hitseb, axis=axis),
            'hitsee': ak.sum(self.hitsee, axis=axis),
            'nevents': ak.sum(self.nevents, axis=axis),
            'nlumis': ak.sum(self.nlumis, axis=axis),
            'fill': ak.min(self.fill, axis=axis),
            'reclumi': ak.sum(self.reclumi, axis=axis),
            'delivlumi': ak.sum(self.delivlumi, axis=axis)
        }
        return ak.zip(
            data,
            with_name='Info'
        )

    def add(self, other):
        """Add a new info to existing one, not in place"""
        data = {
            'maxmiseb': other.maxmiseb,
            'maxmisee': other.maxmisee,
            'minmiseb': other.minmiseb,
            'minmisee': other.minmisee,
            'nmis': other.nmis,
            'hitseb': self.hitseb + other.hitseb,
            'hitsee': self.hitsee + other.hitsee,
            'nevents': self.nevents + other.nevents,
            'nlumis': self.nlumis + other.nlumis,
            'fill': ak.min(other.fill),
            'reclumi': self.reclumi + other.reclumi,
            'delivlumi': self.delivlumi + other.delivlumi
        }
        return ak.zip(
            data,
            with_name='Info'
        )
    
__all__ = ['Info']
