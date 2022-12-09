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
    `nhits` are also saved. `sumet_v` is an array that holds the sum Et values after
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
        # transpose only the last two dimensions to reconstruct the proper sumet_v array
        ax = [-2, -1]
        ax.extend([d for d in range(np.zeros_like(fields).ndim-2)])
        return ak.Array(np.transpose(fields, axes=ax))/self.sumet-1
    
    @property
    def sumet_err(self):
        """
        The statistical uncertainty on the sum Et value
        """
        return np.sqrt(self.sumet2/self.nhits - self.sumet/self.nhits*self.sumet/self.nhits,
                       where=self.nhits>0, out=ak.to_np(ak.zeros_like(self.sumet2)))

    @property
    def sumlc_err(self):
        """The statistical uncertainty on the sum of the LC value"""
        return np.sqrt(self.sumlc2/self.nhits - self.sumlc/self.nhits*self.sumlc/self.nhits,
                       where=self.nhits>0, out=ak.to_np(ak.zeros_like(self.sumlc2)))

    def sum(self, axis: int=0, with_name: str='RecHit'):
        """
        Sum an array of RecHits

        :param axis: axis along which the sum is performed.
        :param with_name: summer RecHit collection name.
        """
        data = { f : ak.sum(self[f], axis=axis) for f in dir(self) if 'sumet_m' in f or 'sumet_p' in f }
        data.update({ "id" : ak.min(self.id, axis=axis),
                      "nhits" : ak.sum(self.nhits, axis=axis),
                      "sumet" : ak.sum(self.sumet, axis=axis),
                      "sumet2" : ak.sum(self.sumet2, axis=axis),
                      "sumlc" : ak.sum(self.sumlc, axis=axis),
                      "sumlc2" : ak.sum(self.sumlc2, axis=axis),
                      "status" : ak.max(self.status, axis=axis) })
        return ak.zip(
            data,
            with_name=with_name
        )

    def add(self, other, with_name: str='RecHit'):
        """
        Add two arrays of RecHits, not inplace.

        :param other: another RecHit.
        """
        data = { f : self[f]+other[f] for f in dir(self) if 'sumet_m' in f or 'sumet_p' in f }
        data.update({ "id": other.id,
                      "nhits": self.nhits + other.nhits,
                      "sumet": self.sumet + other.sumet,
                      "sumet2": self.sumet2 + other.sumet2,
                      "sumlc": self.sumlc + other.sumlc,
                      "sumlc2": self.sumlc2 + other.sumlc2,
                      "status": other.status })
        return ak.zip(
            data,
            with_name=with_name
        )
    
@ak.mixin_class(behavior)
class RecHitEB(RecHit):
    """
    A EcalPhiSym RecHit from the Ecal barrel. This class extends RecHit by adding
    the ieta and iphi methods.
    """

    def sum(self, axis: int=0):
        """
        Sum an array of RecHitEB
        
        :param axis: axis along which the sum is performed.
        """

        return super().sum(axis=axis, with_name='RecHitEB')

    def add(self, other):
        return super().add(other, with_name='RecHitEB')
    
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

    def sum(self, axis: int=0):
        """
        Sum an array of RecHitEE
        
        :param axis: axis along which the sum is performed.
        """

        return super().sum(axis=axis, with_name='RecHitEE')

    def add(self, other):
        return super().add(other, with_name='RecHitEE')
    
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
