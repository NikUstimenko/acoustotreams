"""Scalar basis sets and core array functionalities."""

import abc
from collections import namedtuple

import numpy as np

import treams._operators as op
import acoustotreams._operatorsacoustics as opa
import treams.lattice as la
from treams import util
from treams._lattice import Lattice, WaveVector
from acoustotreams._materialacoustics import AcousticMaterial

class ScalarBasisSet(util.OrderedSet, metaclass=abc.ABCMeta):
    """Parent class of scalar basis sets.

    It is the parent class for all basis sets used. They are expected to be an ordered
    sequence of the modes included in an expansion. Basis sets are expected to
    be immutable.
    """

    _names = ()
    """Names of the relevant parameters"""

    def __repr__(self):
        """String representation.

        Automatically generated when the attribute ``_names`` is defined.

        Returns:
            str
        """
        string = ",\n    ".join(f"{name}={i}" for name, i in zip(self._names, self[()]))
        return f"{self.__class__.__name__}(\n    {string},\n)"

    @classmethod
    @abc.abstractmethod
    def default(cls, *args, **kwargs):
        """Construct a default basis from parameters.

        Construct a basis set in a default order by giving few parameters.
        """
        raise NotImplementedError


class ScalarSphericalWaveBasis(ScalarBasisSet):

    r"""Basis of scalar spherical waves.

    Functions of the spherical wave basis are defined by their degree ``l`` 
    and order ``m``. If the basis is defined with respect to a single 
    expansion center, it is referred to as "global". If it contains multiple 
    expansion centers, it is referred to as "local". In a local basis, an
    additional position index ``pidx`` is used to link the modes to one of the
    specified ``positions``.

    Spherical waves can be separated into regular and singular modes, which
    describe incident and scattered fields, respectively. Therefore, 
    the basis modes refer to one of the functions, :func:`~acoustotreams.special.ssw_rPsi` 
    or :func:`~acoustotreams.special.ssw_Psi`, respectively.

    Args:
        modes (array_like): A tuple containing a list for each of ``l`` and ``m``,
            or ``pidx``, ``l`` and ``m``.
        positions (array_like, optional): The positions of the expansion centers 
            for the specified modes (m). Defaults to ``[[0, 0, 0]]``.

    Attributes:
        pidx (array_like): Integer referring to a row in :attr:`positions`.
        l (array_like): Degree as an integer :math:`l \geq 0`
        m (array_like): Order as an integer :math:`m \leq |l|`.
    """

    _names = ("pidx", "l", "m")

    def __init__(self, modes, positions=None):
        """Initalization."""
        tmp = []
        if isinstance(modes, np.ndarray):
            modes = modes.tolist()
        for m in modes:
            if m not in tmp:
                tmp.append(m)
        modes = tmp
        if len(modes) == 0:
            pidx = []
            l = []  # noqa: E741
            m = []
        elif len(modes[0]) == 2:
            l, m = (*zip(*modes),)
            pidx = np.zeros_like(l)
        elif len(modes[0]) == 3:
            pidx, l, m = (*zip(*modes),)
        else:
            raise ValueError("invalid shape of modes")

        if positions is None:
            positions = np.zeros((1, 3))
        positions = np.array(positions, float)
        if positions.ndim == 1:
            positions = positions[None, :]
        if positions.ndim != 2 or positions.shape[1] != 3:
            raise ValueError(f"invalid shape of positions {positions.shape}")

        self.pidx, self.l, self.m = [
            np.array(i, int) for i in (pidx, l, m)
        ]
        for i, j in ((self.pidx, pidx), (self.l, l), (self.m, m)):
            i.flags.writeable = False
            if np.any(i != j):
                raise ValueError("parameters must be integer")
        if np.any(self.l < 0):
            raise ValueError("'l' must be a non-negative integer")
        if np.any(self.l < np.abs(self.m)):
            raise ValueError("'|m|' cannot be larger than 'l'")
        if np.any(self.pidx >= len(positions)):
            raise ValueError("undefined position is indexed")

        self._positions = positions
        self._positions.flags.writeable = False
        self.lattice = self.kpar = None 

    def __len__(self):
        """Number of modes."""
        return len(self.l)     

    @property
    def positions(self):
        """Positions of the modes' expansion centers.

        The positions are an immutable (N, 3)-array. Each row corresponds to a point in
        the three-dimensional Cartesian space.
        """
        return self._positions

    def __repr__(self):
        """String representation."""
        positions = "positions=" + str(self.positions).replace("\n", ",")
        return f"{super().__repr__()[:-1]}    {positions},\n)"

    @property
    def isglobal(self):
        """Basis is defined with respect to a single (global) expansion center.

        Returns:
            bool
        """
        return len(self) == 0 or np.all(self.pidx == self.pidx[0])

    def __getattr__(self, key):
        dct = {"l": "l", "m": "m", "p": "pidx"}
        try:
            return tuple(getattr(self, dct[k]) for k in key)
        except KeyError:
            raise AttributeError(
                f"'{type(self).__name__}' object has no attribute '{key}'"
            ) from None

    def __getitem__(self, idx):
        """Get a subset of the basis.

        This function supports indexing the basis set using an integer, 
        a slice, an ellipsis, an empty tuple, or a sequence of integers 
        or booleans. Indexing with an integer or an empty tuple returns 
        a tuple; all other index types return a new basis set.

        Alternatively, the string "plm", or "lm" can be used to select 
        only a corresponding subset of the attrubutes :attr:`pidx`, :attr:`l`, 
        and :attr:`m`.
        """
        res = self.pidx[idx], self.l[idx], self.m[idx]
        if isinstance(idx, (int, np.integer)) or (
            isinstance(idx, tuple) and len(idx) == 0
        ):
            return res
        return type(self)(zip(*res), self.positions)

    def __eq__(self, other):
        """Compare basis sets.

        Two basis sets are considered equal if they contain the same modes 
        in the same order and share identical expansion centers :attr:`positions`.
        """
        try:
            return self is other or (
                np.array_equal(self.pidx, other.pidx)
                and np.array_equal(self.l, other.l)
                and np.array_equal(self.m, other.m)
                and np.array_equal(self.positions, other.positions)
            )
        except AttributeError:
            return False

    @classmethod
    def default(cls, lmax, nmax=1, positions=None):
        """Default basis for a given maximum multipolar degree.

        By default, modes are grouped into separate blocks for each position 
        index in ascending order. Within each block, modes are sorted by 
        degree :math:`l` (lowest first) and then by order :math:m from 
        :math:`m = -l` to :math:`m = l`.

        Example:
            >>> acoustotreams.ScalarSphericalWaveBasis.default(2)
            ScalarSphericalWaveBasis(
                pidx=[0 0 0 0 0 0 0 0 0],
                l=[0 1 1 1 2 2 2 2 2],
                m=[ 0 -1  0  1 -2 -1  0  1  2],
                positions=[[0. 0. 0.]],
            )
            >>> acoustotreams.ScalarSphericalWaveBasis.default(0, 2, [[0, 0, 1.], [0, 0, -1.]])
            ScalarSphericalWaveBasis(
                pidx=[0 1],
                l=[0 0],
                m=[0 0],
                positions=[[ 0.  0.  1.], [ 0.  0. -1.]],
            )

        Args:
            lmax (int): Maximum multipolar degree.
            nmax (int, optional): Number of positions. Defaults to 1.
            positions (array_like of float, optional): Positions of the expansion centers (m). Defaults to [[0, 0, 0]].
        """
        modes = [
            [n, l, m]
            for n in range(0, nmax)
            for l in range(0, lmax + 1)  # noqa: E741
            for m in range(-l, l + 1)
        ]
        return cls(modes, positions=positions)

    @staticmethod
    def defaultlmax(dim, nmax=1):
        """Calculate the maximum degree for a given number of modes.

        Return the estimated maximum value of `l` given the dimension of the T-matrix. 
        This is the inverse of :meth:defaultdim. A value of zero is returned for empty T-matrices.

        Example:
            >>> acoustotreams.ScalarSphericalWaveBasis.defaultlmax(len(acoustotreams.ScalarSphericalWaveBasis.default(3)))
            3

        Args:
            dim (int): Dimension of the T-matrix, respectively number of modes.
            nmax (int, optional): Number of scatterers. Defaults to 1.

        Returns:
            int
        """
        res = np.sqrt(dim / nmax) - 1
        res_int = int(np.rint(res))
        if np.abs(res - res_int) > 1e-8 * np.maximum(np.abs(res), np.abs(res_int)):
            raise ValueError("cannot estimate the default lmax")
        return res_int

    @staticmethod
    def defaultdim(lmax, nmax=1):
        """Calculate the number of modes for a given mulipolar degree.

        Return the dimesion of the T-matrix given the maximum value of `l`. 
        This is the inverse of :meth:`defaultlmax`. A value of lmax=-1 is allowed.

        Args:
            lmax (int): Maximum multipolar degree
            nmax (int, optional): Number of scatterers. Defaults to `1`

        Returns:
            int
        """
        # lmax=-1 is allowed and won't give an error
        if lmax < -1 or nmax < 0:
            raise ValueError("maximal degree must be non-negative")
        return (lmax + 1) * (lmax + 1) * nmax

    @classmethod
    def _from_iterable(cls, it, positions=None):
        if isinstance(cls, ScalarSphericalWaveBasis):
            positions = cls.positions if positions is None else positions
            cls = type(cls)
        obj = cls(it, positions=positions)
        return obj    
    

class ScalarCylindricalWaveBasis(ScalarBasisSet):
    r"""Basis of scalar cylindrical waves.

    Functions of the cylindrical wave basis are defined by the Z-components of the wave
    vector ``kz`` and the order  ``m``. If the basis is defined with respect to a single 
    expansion center, it is referred to as "global". If it contains multiple expansion
    centers, it is referred to as "local". In a local basis, an additional position index 
    ``pidx`` is used to link the modes to one of the specified ``positions``.

    Cylindrical waves can be separated into regular and singular modes, which describe 
    incident and scattered fields, respectively. Thus, the basis modes refer to one of 
    the functions :func:`~acoustotreams.special.scw_rPsi` or :func:`~acoustotreams.special.scw_Psi`,
    respectively.

    Args:
        modes (array_like): A tuple containing a list for each of ``kz`` and ``m``,
            or ``pidx``, ``kz`` and ``m``.
        positions (array_like, optional): Positions of the expansion centers for the specified
            modes (m). Defaults to ``[[0, 0, 0]]``.

    Attributes:
        pidx (array_like): Integer referring to a row in :attr:`positions`.
        kz (array_like): Real-valued Z-component of the wave vector (rad/m).
        m (array_like): Integer order.
    """

    _names = ("pidx", "kz", "m")

    def __init__(self, modes, positions=None):
        """Initalization."""
        tmp = []
        if isinstance(modes, np.ndarray):
            modes = modes.tolist()
        for m in modes:
            if m not in tmp:
                tmp.append(m)
        modes = tmp
        if len(modes) == 0:
            pidx = []
            kz = []
            m = []
        elif len(modes[0]) == 2:
            kz, m = (*zip(*modes),)
            pidx = np.zeros_like(m)
        elif len(modes[0]) == 3:
            pidx, kz, m = (*zip(*modes),)
        else:
            raise ValueError("invalid shape of modes")

        if positions is None:
            positions = np.zeros((1, 3))
        positions = np.array(positions, float)
        if positions.ndim == 1:
            positions = positions[None, :]
        if positions.ndim != 2 or positions.shape[1] != 3:
            raise ValueError(f"invalid shape of positions {positions.shape}")

        self.pidx, self.m = [np.array(i, int) for i in (pidx, m)]
        self.kz = np.array(kz, float)
        self.kz.flags.writeable = False
        for i, j in ((self.pidx, pidx), (self.m, m)):
            i.flags.writeable = False
            if np.any(i != j):
                raise ValueError("parameters must be integer")
        if np.any(self.pidx >= len(positions)):
            raise ValueError("undefined position is indexed")

        self._positions = positions
        self._positions.flags.writeable = False

        self.lattice = self.kpar = None
        if len(self.kz) > 0 and np.all(self.kz == self.kz[0]):
            self.kpar = WaveVector(self.kz[0])
    
    def __len__(self):
        """Number of modes."""
        return len(self.m) 

    @property
    def positions(self):
        """Positions of the modes' expansion centers.

        The positions are an immutable (N, 3)-array. Each row corresponds to a point in
        the three-dimensional Cartesian space.
        """
        return self._positions

    def __repr__(self):
        """String representation."""
        positions = "positions=" + str(self.positions).replace("\n", ",")
        return f"{super().__repr__()[:-1]}    {positions},\n)"

    @property
    def isglobal(self):
        """Basis is defined with respect to a single (global) expansion center.

        Returns:
            bool
        """
        return len(self) == 0 or np.all(self.pidx == self.pidx[0])

    def __getattr__(self, key):
        dct = {"z": "kz", "m": "m", "p": "pidx"}
        try:
            return tuple(getattr(self, dct[k]) for k in key)
        except KeyError:
            raise AttributeError(
                f"'{type(self).__name__}' object has no attribute '{key}'"
            ) from None

    def __getitem__(self, idx):
        """Get a subset of the basis.

        This function supports indexing the basis set using an integer, a slice, 
        an ellipsis, an empty tuple, or a sequence of integers or booleans. Indexing 
        with an integer or an empty tuple returns a tuple; all other index types 
        return a new basis set.

        Alternatively, the string "pkzm" or "kzm" can be used to select a 
        corresponding subset of the attributes :attr:`pidx`, :attr:`kz`, and :attr:`m`.
        """
        res = self.pidx[idx], self.kz[idx], self.m[idx]
        if isinstance(idx, (int, np.integer)) or (isinstance(idx, tuple) and idx == ()):
            return res
        return type(self)(zip(*res), self.positions)

    def __eq__(self, other):
        """Compare basis sets.

        Two basis sets are considered equal if they contain the same modes in the same 
        order and share identical expansion centers :attr:`positions`.
        """
        try:
            return self is other or (
                np.array_equal(self.pidx, other.pidx)
                and np.array_equal(self.kz, other.kz)
                and np.array_equal(self.m, other.m)
                and np.array_equal(self.positions, other.positions)
            )
        except AttributeError:
            return False

    @classmethod
    def default(cls, kzs, mmax, nmax=1, positions=None):
        """Default basis for given Z-components of the wave vector and a given maximum order.

        By default, modes are grouped into separate blocks for each position 
        index in ascending order. Within each block, modes are sorted by 
        the Z-component of the wave vector :math:`k_z` and then by order :math:m 
        in asceding order.

        Example:
            >>> acoustotreams.ScalarCylindricalWaveBasis.default([-0.5, 0.5], 1)
            ScalarCylindricalWaveBasis(
                pidx=[0 0 0 0 0 0],
                kz=[-0.5 -0.5 -0.5 0.5  0.5  0.5],
                m=[-1  0  1 -1  0  1],
                positions=[[0. 0. 0.]],
            )
            >>> acoustotreams.ScalarCylindricalWaveBasis.default([0], 1, 2, [[1., 0, 0], [-1., 0, 0]])
            ScalarCylindricalWaveBasis(
                pidx=[0 0 0 1 1 1],
                kz=[0. 0. 0. 0. 0. 0.],
                m=[-1  0  1 -1  0  1],
                positions=[[ 1.  0.  0.], [-1.  0.  0.]],
            )

        Args:
            kzs (array_like of float): Z-components of the wave vector (rad/m).
            mmax (int): Maximum value of order.
            nmax (int, optional): Number of positions. Defaults to 1.
            positions (array_like of float, optional): Positions of the expansion centers (m). Defaults to [[0, 0, 0]].
        """
        kzs = np.atleast_1d(kzs)
        if kzs.ndim > 1:
            raise ValueError(f"kzs has dimension larger than one: '{kzs.ndim}'")
        modes = [
            [n, kz, m]
            for n in range(nmax)
            for kz in kzs
            for m in range(-mmax, mmax + 1)
        ]
        return cls(modes, positions=positions)

    @classmethod
    def diffr_orders(cls, kz, mmax, lattice, bmax, nmax=1, positions=None):
        """Creates a basis set for a system periodic in the z-direction.

        Example:
           >>> acoustotreams.ScalarCylindricalWaveBasis.diffr_orders(0.1, 1, 2 * np.pi, 1)
           ScalarCylindricalWaveBasis(
               pidx=[0 0 0 0 0 0 0 0 0],
               kz=[-0.9 -0.9 -0.9  0.1  0.1  0.1  1.1  1.1  1.1],
               m=[-1  0  1 -1  0  1 -1  0  1],
               positions=[[0. 0. 0.]],
           )

        Args:
            kz (float): Z-component of the wave vector. Ideally, it is defined within the 
                first Brillouin zone (use :func:`treams.misc.firstbrillouin1d`).
            mmax (int): Maximum value of the order.
            lattice (:class:`acoustotreams.Lattice` or float): Lattice definition or pitch.
            bmax (float): Maximum length of reciprocal-lattice vectors that defines 
                a maximum momentum transfer from the given value of `kz`.
            nmax (int, optional): Number of expansion centers.
            positions (array_like of float, optional): Positions of the expansion centers (m). Defaults to [[0, 0, 0]].
        """
        lattice = Lattice(lattice)
        lattice_z = Lattice(lattice, "z")
        nkz = np.floor(np.abs(bmax / lattice_z.reciprocal))
        kzs = kz + np.arange(-nkz, nkz + 1) * lattice_z.reciprocal
        res = cls.default(kzs, mmax, nmax, positions=positions)
        res.lattice = lattice
        res.kpar = WaveVector(kz)
        return res

    @classmethod
    def _from_iterable(cls, it, positions=None):
        if isinstance(cls, ScalarCylindricalWaveBasis):
            positions = cls.positions if positions is None else positions
            lattice = cls.lattice
            kpar = cls.kpar
            cls = type(cls)
        else:
            lattice = kpar = None
        obj = cls(it, positions=positions)
        obj.lattice = lattice
        obj.kpar = kpar
        return obj

    @staticmethod
    def defaultmmax(dim, nkz=1, nmax=1):
        """Calculate the maximum order for a given number of modes.

        Return the estimated maximum value of `m` given the dimension of the T-matrix.
        This is the inverse of :meth:`defaultdim`. A value of zero is returned for 
        empty T-matrices.

        Example:
            >>> acoustotreams.ScalarCylindricalWaveBasis.defaultmmax(len(acoustotreams.ScalarCylindricalWaveBasis.default([0], 2)), 1)
            2

        Args:
            dim (int): Dimension of the T-matrix, respectively number of modes.
            nkz (int, optional): Number of Z-components of the wave vector. Defaults to 1.
            nmax (int, optional): Number of scatterers. Defaults to 1.

        Returns:
            int
        """
        if dim % (nkz * nmax) != 0:
            raise ValueError("cannot estimate the default mmax")
        dim = dim // (nkz * nmax)
        if (dim - 1) % 2 != 0:
            raise ValueError("cannot estimate the default mmax")
        return (dim - 1) // 2

    @staticmethod
    def defaultdim(nkz, mmax, nmax=1):
        """Calculate the number of modes for a given mulipolar order.

        Return the dimension of the T-matrix given the maximum value of 
        `m` and number of `k_z` components. A value of zero is allowed.

        Args:
            nkz (int): Number of Z-components of the wave vector.
            mmax (int): Maximum value of the order.
            nmax (int, optional): Number of scatterers. Defaults to 1.

        Returns:
            int
        """
        if nkz < 0 or mmax < 0:
            raise ValueError("maximal order must be positive")
        return (2 * mmax + 1) * nkz * nmax
    
    
class ScalarPlaneWaveBasis(ScalarBasisSet):
    """Scalar plane wave basis parent class."""

    isglobal = True


class ScalarPlaneWaveBasisByUnitVector(ScalarPlaneWaveBasis):
    """Scalar plane wave basis.

    A plane-wave basis is defined by a collection of wave vectors specified by their
    Cartesian components ``qx``, ``qy``, and ``qz`` normalized to
    :math:`q_x^2 + q_y^2 + q_z^2 = 1`.

    Plane waves can refer to:func:`~acoustotreams.special.spw_Psi`.

    Args:
        modes (array_like): A tuple containing a list for each of ``qx``, ``qy``,
            and ``qz``.

    Attributes:
        qx (array_like): X-component of the normalized wave vector.
        qy (array_like): Y-component of the normalized wave vector.
        qz (array_like): Z-component of the normalized wave vector.
    """

    _names = ("qx", "qy", "qz")
    """A scalar plane-wave basis is always global."""

    def __init__(self, modes):
        """Initialization."""
        tmp = []
        if isinstance(modes, np.ndarray):
            modes = modes.tolist()
        for m in modes:
            if m not in tmp:
                tmp.append(m)
        modes = tmp
        if len(modes) == 0:
            qx = []
            qy = []
            qz = []
        elif len(modes[0]) == 3:
            qx, qy, qz = (*zip(*modes),)
        else:
            raise ValueError("invalid shape of modes")

        qx, qy, qz = (np.array(i) for i in (qx, qy, qz))
        norm = np.emath.sqrt(qx * qx + qy * qy + qz * qz)
        norm[np.abs(norm - 1) < 1e-14] = 1
        qx, qy, qz = (np.true_divide(i, norm) for i in (qx, qy, qz))
        for i in (qx, qy, qz):
            i.flags.writeable = False
            if i.ndim > 1:
                raise ValueError("invalid shape of parameters")
        self.qx, self.qy, self.qz = qx, qy, qz
        self.lattice = self.kpar = None

    def __len__(self):
        """Number of modes."""
        return len(self.qx) 

    def __getattr__(self, key):
        dct = {"x": "qx", "y": "qy", "z": "qz"}
        try:
            return tuple(getattr(self, dct[k]) for k in key)
        except KeyError:
            raise AttributeError(
                f"'{type(self).__name__}' object has no attribute '{key}'"
            ) from None

    def __getitem__(self, idx):
        """Get a subset of the basis.

        This function supports indexing the basis set using an integer, a slice, an ellipsis, 
        an empty tuple, or a sequence of integers or booleans. Indexing with an integer or 
        an empty tuple returns a tuple; all other index types return a new basis set.

        Alternatively, the string "xyz", "xy" or "z" can be used to select only a
        corresponding subset of the attributes :attr:`qx`, :attr:`qy`, and :attr:`qz`.
        """
        res = self.qx[idx], self.qy[idx], self.qz[idx]
        if isinstance(idx, (int, np.integer)) or (isinstance(idx, tuple) and idx == ()):
            return res
        return type(self)(zip(*res))

    @classmethod
    def default(cls, kvecs):
        """Default basis from the given wave vectors.

        Example:
            >>> acoustotreams.ScalarPlaneWaveBasisByUnitVector.default([[0, 0, 5], [0, 3, 4]])
            ScalarPlaneWaveBasisByUnitVector(
                qx=[0. 0.],
                qy=[0.  0.6],
                qz=[1.  0.8],
            )

        Args:
            kvecs (array_like): Wave vectors in Cartesian coordinates (rad/m).
        """
        kvecs = np.atleast_2d(kvecs)
        modes = np.empty((kvecs.shape[0], 3), kvecs.dtype)                           
        modes = kvecs
        return cls(modes)

    @classmethod
    def _from_iterable(cls, it):
        if isinstance(cls, ScalarPlaneWaveBasisByUnitVector):
            lattice = cls.lattice
            kpar = cls.kpar
            cls = type(cls)
        else:
            lattice = kpar = None
        obj = cls(it)
        obj.lattice = lattice
        obj.kpar = kpar
        return obj

    def __eq__(self, other):
        """Compares basis sets.

        Two basis sets are considered equal if they contain the same modes in the same order.
        """
        try:
            return self is other or (
                np.array_equal(self.qx, other.qx)
                and np.array_equal(self.qy, other.qy)
                and np.array_equal(self.qz, other.qz)
            )
        except AttributeError:
            return False

    def bycomp(self, k0, alignment="xy", material=AcousticMaterial()):
        """Creates an instance of :class:`ScalarPlaneWaveBasisByComp`.

        The plane-wave basis is changed to a partial basis, where only two (real-valued)
        wave-vector components are specified. A third component is then inferred from the
        dispersion relation, which depends on the wavenumber, material and
        ``modetype``, which determines the sign of the third component.

        Args:
            alignment (str, optional): Wave-vector components included in the
                partial basis. Defaults to "xy", other permitted values are "yz" and
                "zx".
            k0 (float, optional): Angular wavenumber in the air (rad/m). 
                Verifies that the current basis satisfies the dispersion relation.
            material (:class:`acoustotreams.AcousticMaterial` or tuple): Material definition. 
                Defaults to the air.
        """
        ks = material.ks(k0)
        if alignment in ("xy", "yz", "zx"):
            kpars = [ks * getattr(self, "q" + s) for s in alignment]
        else:
            raise ValueError(f"invalid alignment '{alignment}'")
        obj = ScalarPlaneWaveBasisByComp(zip(*kpars), alignment)
        obj.lattice = self.lattice
        obj.kpar = self.kpar
        return obj

    def kvecs(self, k0, material=AcousticMaterial(), modetype=None):
        """Wave vectors.

        Args:
            k0 (float): Angular wavenumber in the air (rad/m).
            material (:class:`acoustotreams.AcousticMaterial` or tuple, optional): Material 
                definition. Defaults to the air.
            modetype (optional): Currently unused for this class.
        """
        # TODO: check kz depending on modetype (alignment?)
        ks = AcousticMaterial(material).ks(k0)
        return ks * self.qx, ks * self.qy, ks * self.qz

    def permute(self, n=1):
        n = n % 3
        lattice = None if self.lattice is None else self.lattice.permute(n)
        kpar = None if self.kpar is None else self.kpar.permute(n)
        qx, qy, qz = self.qx, self.qy, self.qz
        for _ in range(n):
            qx, qy, qz = qz, qx, qy
        obj = type(self)(zip(*(qx, qy, qz)))
        obj.lattice = lattice
        obj.kpar = kpar
        return obj
    

class ScalarPlaneWaveBasisByComp(ScalarPlaneWaveBasis):
    """Partial scalar plane wave basis.

    A partial plane-wave basis is defined by two wave-vector components out of all three
    Cartesian wave-vector components ``kx``, ``ky``, and ``kz``. 
    Which two components are given is specified in the :attr:`alignment`. This
    basis is mostly used for stratified media that are periodic or uniform in the two
    alignment directions, such that the given wave-vector components correspond to the
    diffraction orders.

    Scalar plane waves refer to :func:`~acoustotreams.special.spw_Psi`.

    Args:
        modes (array_like): A tuple containing a list for each of ``k1`` and ``k2``.
        alignment (str, optional): Wave-vector components included in the
                partial basis. Defaults to "xy", other possible values are "yz" and "zx".

    Attributes:
        alignment (str): Alignment of the partial basis.
    """

    def __init__(self, modes, alignment="xy"):
        """Initialization."""
        if isinstance(modes, np.ndarray):
            modes = modes.tolist()
        tmp = []
        for m in modes:
            if m not in tmp:
                tmp.append(m)
        modes = tmp
        if len(modes) == 0:
            kx = []
            ky = []
        elif len(modes[0]) == 2:
            kx, ky = (*zip(*modes),)
        else:
            raise ValueError("invalid shape of modes")

        self._kx, self._ky = [np.real(i) for i in (kx, ky)]
        for i, j in [(self._kx, kx), (self._ky, ky)]:
            i.flags.writeable = False
            if i.ndim > 1:
                raise ValueError("invalid shape of parameters")
            if np.any(i != j):
                raise ValueError("invalid value for parameter, must be real")

        if alignment in ("xy", "yz", "zx"):
            self._names = (*("k" + i for i in alignment),)
        else:
            raise ValueError(f"invalid alignment '{alignment}'")

        self.alignment = alignment
        self.lattice = self.kpar = None

    def __len__(self):
        """Number of modes."""
        return len(self.kx)

    def permute(self, n=1):
        n = n % 3
        lattice = None if self.lattice is None else self.lattice.permute(n)
        kpar = None if self.kpar is None else self.kpar.permute(n)
        alignments = {"xy": "yz", "yz": "zx", "zx": "xy"}
        alignment = self.alignment
        for _ in range(n):
            alignment = alignments[alignment]
        obj = self._from_iterable(self, alignment=alignment)
        obj.lattice = lattice
        obj.kpar = kpar
        return obj

    @property
    def kx(self):
        """X-components of the wave vector.

        If the components are not specified, `None` is returned.
        """
        if self.alignment == "xy":
            return self._kx
        if self.alignment == "zx":
            return self._ky
        return None

    @property
    def ky(self):
        """Y-components of the wave vector.

        If the components are not specified, `None` is returned.
        """
        if self.alignment == "xy":
            return self._ky
        if self.alignment == "yz":
            return self._kx
        return None

    @property
    def kz(self):
        """Z-components of the wave vector.

        If the components are not specified, `None` is returned.
        """
        if self.alignment == "yz":
            return self._ky
        if self.alignment == "zx":
            return self._kx
        return None

    def __getattr__(self, key):
        dct = {"x": "kx", "y": "ky", "z": "kz"}
        try:
            return tuple(getattr(self, dct[k]) for k in key)
        except KeyError:
            raise AttributeError(
                f"'{type(self).__name__}' object has no attribute '{key}'"
            ) from None

    def __getitem__(self, idx):
        """Get a subset of the basis.

        This function supports indexing the basis set using an integer, a slice, 
        an ellipsis, an empty tuple, or a sequence of integers or booleans. 
        Indexing with an integer or an empty tuple returns a tuple; 
        all other index types return a new basis set.
        """
        res = self._kx[idx], self._ky[idx]
        if isinstance(idx, (int, np.integer)) or (isinstance(idx, tuple) and idx == ()):
            return res
        return type(self)(zip(*res))

    @classmethod
    def default(cls, kpars, alignment="xy"):
        """Default basis from given wave vectors.

        Example:
            >>> acoustotreams.ScalarPlaneWaveBasisByComp.default([[0, 0], [0, 3]])
            ScalarPlaneWaveBasisByComp(
                kx=[0 0],
                ky=[0 3],
            )

        Args:
            kpars (array_like): Two Cartesian components of the wave vectors (rad/m).
            alignment (str, optional): Wave-vector components included in the
                partial basis. Defaults to "xy", other possible values are "yz" and "zx".

        """
        kpars = np.atleast_2d(kpars)
        modes = kpars
        return cls(modes, alignment=alignment)

    @classmethod
    def diffr_orders(cls, kpar, lattice, bmax):
        """Creates a basis set for a two-dimensional periodic system.

        The reciprocal lattice of the given lattice is used to determine 
        all diffraction orders within the defined maximum radius of the circle 
        in the reciprocal space.

        Example:
            >>> acoustotreams.ScalarPlaneWaveBasisByComp.diffr_orders([0, 0], acoustotreams.Lattice.square(2 * np.pi), 1)
            ScalarPlaneWaveBasisByComp(
                kx=[ 0.  0.  0.  1. -1.],
                ky=[ 0.  1. -1.  0.  0.],
            )

        Args:
            kpar (float): In-plane wave-vector components (rad/m). Ideally, they are defined within the
                first Brillouin zone (use :func:`misc.firstbrillouin2d`).
            lattice (:class:`acoustotreams.Lattice` or float): Lattice definition or pitch.
            bmax (float): Maximum length of reciprocal-lattice vectors that defines a maximum
                momentum transfer with respect to `kpar`.
        """
        lattice = Lattice(lattice)
        if lattice.dim != 2:
            raise ValueError("invalid lattice dimensions")
        latrec = lattice.reciprocal
        kpars = kpar + la.diffr_orders_circle(latrec, bmax) @ latrec
        obj = cls.default(kpars, alignment=lattice.alignment)
        obj.lattice = lattice
        obj.kpar = WaveVector(kpar, alignment=lattice.alignment)
        return obj

    @classmethod
    def _from_iterable(cls, it, alignment="xy"):
        if isinstance(cls, ScalarPlaneWaveBasisByComp):
            alignment = cls.alignment if alignment is None else alignment
            lattice = cls.lattice
            kpar = cls.kpar
            cls = type(cls)
        else:
            lattice = kpar = None
        obj = cls(it, alignment)
        obj.lattice = lattice
        obj.kpar = kpar
        return obj

    def __eq__(self, other):
        """Compares basis sets.

        Two basis sets are considered equal if they contain the same modes in the same order
        """
        try:
            skx, sky, skz = self.kx, self.ky, self.kz
            okx, oky, okz = other.kx, other.ky, other.kz
            return self is other or (
                (np.array_equal(skx, okx) or (skx is None and okx is None))
                and (np.array_equal(sky, oky) or (sky is None and oky is None))
                and (np.array_equal(skz, okz) or (skz is None and okz is None))
            )
        except AttributeError:
            return False

    def byunitvector(self, k0, material=AcousticMaterial(), modetype="up"):
        """Creates an instance of complete basis :class:`ScalarPlaneWaveBasis`.

        A plane-wave basis is considered complete when all three Cartesian components
        are defined for each mode. A third Cartesian component is calculated from
        the specified wavenumber, material, and modetype. The modetype "up" 
        (resp. "down") corresponds to waves propagating in the positive (resp. negative)
        direction with respect to the Cartesian axis that is orthogonal to the
        alignment plane.

        Args:
            k0 (float): Angular wavenumber in the air (rad/m).
            material (:class:`~acoustotreams.AcousticMaterial` or tuple, optional): Material 
                definition. Defaults to the air.
            modetype (str, optional): Propagation direction. Defaults to "up".
        """
        if modetype not in ("up", "down"):
            raise ValueError("modetype not recognized")
        material = AcousticMaterial(material)
        kx = self._kx
        ky = self._ky
        kz = material.kzs(k0, kx, ky) * (2 * (modetype == "up") - 1)
        if self.alignment == "yz":
            kx, ky, kz = kz, kx, ky
        elif self.alignment == "zy":
            kx, ky, kz = ky, kz, kx
        obj = ScalarPlaneWaveBasisByUnitVector(zip(kx, ky, kz))
        obj.lattice = self.lattice
        obj.kpar = self.kpar
        return obj

    def kvecs(self, k0, material=AcousticMaterial(), modetype="up"):
        """Wave vectors.

        Args:
            k0 (float): Angular wavenumber in the air (rad/m).
            material (:class:`~acoustotreams.AcousticMaterial` or tuple, optional): Material
                definition. Defaults to the air.
            modetype (str, optional): Propagation direction. Defaults to "up".
        """
        if modetype not in ("up", "down"):
            raise ValueError("modetype not recognized")
        material = AcousticMaterial(material)
        kx = self._kx
        ky = self._ky
        kz = material.kzs(k0, kx, ky) * (2 * (modetype == "up") - 1)
        if self.alignment == "yz":
            return kz, kx, ky
        if self.alignment == "zx":
            return ky, kz, kx
        return kx, ky, kz


def _raise_basis_error(*args):
    raise TypeError("'basis' must be ScalarBasisSet")


class ScalarPhysicsDict(util.AnnotationDict):
    """Physics dictionary (for scalar waves).

    Derives from :class:`treams.util.AnnotationDict`. This dictionary has additionally
    several properties defined.

    Attributes:
        basis (:class:`ScalarBasisSet`): Basis set.
        k0 (float): Angular wavenumber in the air (rad/m).
        kpar (list): In-plane wave-vector components (rad/m). Usually, this is a list of length
            3 with its items corresponding to the Cartesian coordinates. Unspecified items
            aree set to `nan`.
        lattice (:class:`~acoustotreams.Lattice`): Lattice definition.
        material (:class:`~acoustotreams.AcousticMaterial`): Material definition.
        modetype (str): Mode type. For spherical and cylindrical waves, it can be
            "incident" and "scattered"; for partial plane waves, "up" or "down".
    """

    properties = {
        "basis": (
            ":class:`ScalarBasisSet`.",
            lambda x: isinstance(x, ScalarBasisSet),
            _raise_basis_error,
        ),
        "k0": ("Angular wavenumber.", lambda x: isinstance(x, float), float),
        "kpar": (
            "Wave-vector components tangential to the lattice.",
            lambda x: isinstance(x, WaveVector),
            WaveVector,
        ),
        "lattice": (
            ":class:`~acoustotreams.Lattice`.",
            lambda x: isinstance(x, Lattice),
            Lattice,
        ),
        "material": (
            ":class:`AcousticMaterial`.",
            lambda x: isinstance(x, AcousticMaterial),
            AcousticMaterial,
        ),
        "modetype": ("Mode type.", lambda x: isinstance(x, str), str),
    }
    """Special properties tracked by the ScalarPhysicsDict."""

    def __setitem__(self, key, val):
        """Sets item specified by key to the defined value.

        When overwriting an existing key an :class:`AnnotationWarning` is emitted.
        Avoid the warning by explicitly deleting the key first. The special attributes
        are cast to their corresponding types.

        Args:
            key (hashable): Key
            val : Value

        Warns:
            AnnotationWarning
        """
        if key not in self.properties:
            raise AttributeError(f"invalid key '{key}'")
        _, testfunc, castfunc = self.properties[key]
        if not testfunc(val):
            val = castfunc(val)
        super().__setitem__(key, val)


class AcousticsArray(util.AnnotatedArray):
    """Acoustics-aware array.

    An acoustics-aware array is a special type of :class`~treams.util.AnnotatedArray`.
    Additionally to keeping track of the annotations, it is enhanced by the ability to
    create suiting linear operators to perform tasks like rotations, translations, or
    expansions into different basis sets and by applying special rules for the
    physical properties of :class:`ScalarPhysicsDict` upon matrix multiplications, see also
    :meth:`__matmul__`.
    """

    _scales = {"basis"}

    pfield = op.OperatorAttribute(opa.PField)
    """Pressure field evaluation matrix, see also :class:`acoustotreams.PField`."""
    vfield = op.OperatorAttribute(opa.VField)
    """Velocity field evaluation matrix, see also :class:`acoustotreams.VField`."""
    pamplitudeff = op.OperatorAttribute(opa.PAmplitudeFF)
    """Far-field amplitude of pressure field evaluation matrix, see also :class:`acoustotreams.PAmplitudeFF`."""
    expand = op.OperatorAttribute(opa.Expand)
    """Expansion matrix, see also :class:`acoustotreams.Expand`."""
    expandlattice = op.OperatorAttribute(opa.ExpandLattice)
    """Lattice expansion matrix, see also :class:`acoustotreams.ExpandLattice`."""
    permute = op.OperatorAttribute(opa.Permute)
    """Permutation matrix, see also :class:`acoustotreams.Permute`."""
    rotate = op.OperatorAttribute(opa.Rotate)
    """Rotation matrix, see also :class:`acoustotreams.Rotate`."""
    translate = op.OperatorAttribute(opa.Translate)
    """Translation matrix, see also :class:`acoustotreams.Translate`."""

    def __init__(self, arr, ann=(), /, **kwargs):
        """Initialization."""
        super().__init__(arr, ann, **kwargs)
        self._check()

    @property
    def ann(self):
        """Array annotations."""
        # necessary to define the setter below
        return super().ann

    def __getattr__(self, key):
        try:
            return super().__getattr__(key)
        except AttributeError as err:
            if key in ScalarPhysicsDict.properties:
                return None
            raise err from None

    def __setattr__(self, key, val):
        if key in ScalarPhysicsDict.properties:
            val = (val,) * self.ndim if not isinstance(val, tuple) else val
            self.ann.as_dict[key] = val
        else:
            super().__setattr__(key, val)

    @ann.setter
    def ann(self, ann):
        """Sets array annotations.

        See also :meth:`treams.util.AnnotatedArray.__setitem__`.
        """
        self._ann = util.AnnotationSequence(*(({},) * self.ndim), mapping=ScalarPhysicsDict)
        self._ann.update(ann)

    def __repr__(self):
        """String representation.

        For a more managable output only the special physics properties are shown
        alongside the array itself.
        """
        repr_arr = "    " + repr(self._array)[5:-1].replace("\n  ", "\n") + ","        
        for key in ScalarPhysicsDict.properties:
            if getattr(self, key) is not None:
                repr_arr += f"\n    {key}={repr(getattr(self, key))},"
        return f"{self.__class__.__name__}(\n{repr_arr}\n)"
    
    def _check(self):
        """Run checks to validate the physical properties.

        The checks are run on the last two dimensions. They include:

            * Dispersion relation checks for :class:`ScalarPlaneWaveBasis` if `basis`, `k0`
              and `material` are defined`
            * All lattices explicitly given and in the basis hints must be compatible.
            * All tangential wave vector compontents must be compatible.
        """
        for a in self.ann[-2:]:
            k0 = a.get("k0")
            material = a.get("material")
            modetype = a.get("modetype")
            basis = a.get("basis", namedtuple("_basis", "lattice kpar")(None, None))

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        """Implements ufunc API.

        In addition to keeping track of the annotations, the special properties of an
        AcousticsArray are also "transparent" in matrix multiplications.
        """
        res = super().__array_ufunc__(ufunc, method, *inputs, **kwargs)
        if (
            ufunc is np.matmul
            and method == "__call__"
            and not isinstance(res, np.generic)
        ):
            axes = kwargs.get(
                "axes", [tuple(range(-min(np.ndim(i), 2), 0)) for i in inputs + (res,)]   #!!!!!!!!!!!
            )
            anns = [
                i.ann[ax] if hasattr(i, "ann") else [{}, {}]
                for i, ax in zip(inputs, axes)
            ]
            for name in ScalarPhysicsDict.properties:
                if name in anns[0][-1] and all(name not in a for a in anns[1]):
                    res.ann[axes[-1][-1]].setdefault(name, anns[0][-1][name])
                if name in anns[1][0] and all(name not in a for a in anns[0]):
                    res.ann[axes[-1][0]].setdefault(name, anns[1][0][name])
        return res