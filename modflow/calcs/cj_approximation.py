import numpy as np


class CooperJacob:

    def __init__(self, q=None, S=None, r=None, t=None, T=None, k=None, b=None, gpm_ft_day=True):
        self.q = q
        self.T = T
        self.S = S
        self.r = r
        self.t = t
        self.k = k
        self.b = b
        self.gpm_ft_day = gpm_ft_day
        if gpm_ft_day and q is not None:
            # convert a gpm Q to a cubic-feet-per-day Q
            self.q = self.q * 192.5


    def cooper_jacob_drawdown(self, q, S, r, t, k=None, b=None, T=None):
        """
        Calculate transient drawdown at an observation well using Cooper-Jacob approximation.
        Can be any consistent units. Unit prompts in documentation are guidance examples.

        Parameters:
        Q (float): Pumping rate at the pumping well (ft³/day)
        T (float): Transmissivity of the aquifer (ft²/day), optional, can provide k and b instead
        S (float): Storage coefficient (dimensionless)
        r (float): Distance from the pumping well to the observation well (ft)
        t (float): Time since pumping started (days)
        k (float): aquifer hydraulic conductivity (ft/day), optional, can provide T instead
        b (float): aquifer thickness (ft), optional, can provide T instead

        Returns:
        float: Drawdown at the observation well (ft)
        """
        q = self.q if q is None else q
        T = self.T if T is None else T
        S = self.S if S is None else S
        r = self.r if r is None else r
        t = self.t if t is None else t
        k = self.k if k is None else k
        b = self.b if b is None else b

        kb = [k, b]
        if T is None:
            for i in kb:
                assert i is not None, f' {i} must be provided!'
            T =  k * b
        qTSrt = [q, T, S, r, t]
        for i in qTSrt:
            assert i is not None, f'{i} must be provided'

        drawdown = ((2.303 * q) / (4 * np.pi * T)) * (np.log((2.25 * T * t) / ((r ** 2) * S)))

        return drawdown

    def get_ds(self, qs, Ss, rs, times, ks=None, bs=None, Ts=None):
        if Ts is None:
            Ts = [k * b for k, b in zip(ks, bs)]
        if self.gpm_ft_day:
            qs = [q * 192.5 for q in qs]
        iterations = zip(qs, Ts, Ss, rs)
        ds_all_inters = []
        for iter in iterations:
            iter_ds = []
            for time in times:
                # get drawdown at time in times for interation
                dt = self.cooper_jacob_drawdown(
                    q=iter[0], T=iter[1], S=iter[2], r=iter[3], t=time)
                iter_ds += [dt]
            ds_all_inters += [iter_ds]
        return ds_all_inters

if __name__ == '__main__':

    # Example usage
    Q_gpm = 1600  # Pumping rate in gpm
    Q_ft3_day = Q_gpm * 192.5  # Convert gpm to ft³/day

    T = 1000  # Transmissivity in ft²/day
    S = 1e-4  # Storage coefficient (dimensionless)
    r = 300  # Distance from the pumping well to the observation well in feet
    t_days = 1  # Time since pumping started in days

    drawdown = cooper_jacob_drawdown(Q_ft3_day, T, S, r, t_days)
    print(f"Drawdown at the observation well: {drawdown:.6f} feet")
