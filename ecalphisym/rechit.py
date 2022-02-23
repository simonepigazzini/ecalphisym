"""EcalPhiSymRecHit with numpy mixins

This provides a class to handle an EcalPhiSym RecHit and
implements high level operations among its data members
"""
import numpy as np
import awkward as ak

behavior = {}

@ak.mixin_class(behavior)
class RecHit:
    """
    A EcalPhiSym RecHit contains information belonging to a single ECAL channel.
    It provide the `sumet` and `sumlc` measurements as well as the squared 
    sums `sumet2` `sumlc2` that are used to compute the statistical uncertainty 
    on the measured values. The channel `status`, `id` and number of recorded hits
    `nhits` are also saved. `sumet_mis` is an array that holds the sum Et values after
    applying a known miscalibration.
    """

    @property
    def sumet_v(self):
        """
        An array containing all the mis-calibratated and nominal sumEts [-m,..,0,..,+m]
        """
        # get miscalibrated sumEt values
        n_mis = len([m for m in dir(self) if 'sumet_m' in m])
        fields = [self['sumet_m'+str(n_mis-i)] for i in range(n_mis)]+[self['sumet']]+[self['sumet_p'+str(i+1)] for i in range(n_mis)]
        return ak.Array(np.transpose(fields, axes=(1, 2, 0)))/self.sumet-1
 
    @property
    def sumet_err(self):
        """
        The statistical uncertainty on the sum Et value
        """
        return np.sqrt(self.sumet2/self.nhits - self.sumet/self.nhits*self.suamet/self.nhits,
                       where=self.nhits>0, out=ak.to_np(ak.zeros_like(self.sumet2)))

    @property
    def sumlc_err(self):
        """The statistical uncertainty on the sum of the LC value"""
        return np.sqrt(self.sumlc2/self.nhits - self.sumlc/self.nhits*self.suamlc/self.nhits,
                       where=self.nhits>0, out=ak.to_np(ak.zeros_like(self.sumlc2)))

    def sum(self, axis: int=0):
        """
        Sum an array of RecHits

        :param axis: axis along which the sum is performed.
        """
        data = { f : ak.sum(self[f], axis=axis) for f in dir(self) if 'sumet_m' in f or 'sumet_p' in f }
        data.update({ "nhits": ak.sum(self.nhits, axis=axis),
                      "sumet": ak.sum(self.sumet, axis=axis),
                      "sumet2": ak.sum(self.sumet2, axis=axis),
                      "sumlc": ak.sum(self.sumlc, axis=axis),
                      "sumlc2": ak.sum(self.sumlc2, axis=axis) })
        return ak.zip(
            data,
            with_name="RecHit"
        )

@ak.mixin_class(behavior)
class RecHitEB(RecHit):
    """
    A EcalPhiSym RecHit from the Ecal barrel. This class extends RecHit by adding
    the ieta and iphi methods.
    """

    def zside(self):
        """
        Barrel eta side
        """
        zside = (self.id & 0x10000)/0x10000
        return 2*zside + ak.full_like(zside, fill_value=-1)
        
    @property
    def ieta(self):
        """
        ieta index of the crystal [-85, 85] (signed)
        """
        return ((self.id >> 9) & 0x7F)*self.zside()

    @property
    def iphi(self):
        """
        iphi index of the crystal [1, 360]
        """
        return self.id & 0x1FF
    
@ak.mixin_class(behavior)
class RecHitEE(RecHit):
    """
    A EcalPhiSym RecHit from the Ecal edcaps. This class extends RecHit by adding
    the ix, iy and iz methods.
    """

    @property
    def iz(self):
        """
        Endcap side
        """
        iz = (self.id & 0x4000)/0x4000
        return 2*iz + ak.full_like(iz, fill_value=-1)
        
    @property
    def ix(self):
        """
        ix index of the crystal [1, 100]
        """
        return (self.id >> 7) & 0x7F

    @property
    def iy(self):
        """
        iy index of the crystal [1, 100]
        """
        return self.id & 0x7F

__all__ = ["RecHit", "RecHitEB", "RecHitEE"]
