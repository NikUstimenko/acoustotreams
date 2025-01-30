import numpy as np
import pytest

import acoustotreams


def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


class TestRotate:
    def test_ssw_invalid(self):
        a = acoustotreams.ScalarSphericalWaveBasis([[1, 0]])
        b = acoustotreams.ScalarCylindricalWaveBasis([[1, -1], [1, 1]])
        with pytest.raises(AttributeError):
            acoustotreams.rotate(1, basis=(b, a))

    def test_ssw(self):
        a = acoustotreams.ScalarSphericalWaveBasis([[1, 0]])
        b = acoustotreams.ScalarSphericalWaveBasis([[1, -1], [1, 1]])
        where = [True, False]
        x = acoustotreams.rotate(1, 2, 3, basis=(a, b), where=where)
        y = acoustotreams.AcousticsArray(
            [[acoustotreams.ssw.rotate(1, 0, 1, -1, 1, 2, 3), 0]], basis=(a, b)
        )
        assert np.all(np.abs(x - y) < 1e-14) and x.ann == y.ann

    def test_scw(self):
        a = acoustotreams.ScalarCylindricalWaveBasis([[0.1, 0]])
        b = acoustotreams.ScalarCylindricalWaveBasis([[0.1, 0], [0.1, 1]])
        where = [True, False]
        x = acoustotreams.rotate(2, basis=(a, b), where=where)
        y = acoustotreams.AcousticsArray(
            [[acoustotreams.scw.rotate(0.1, 0, 0.1, 0, 2), 0]], basis=(a, b)
        )
        assert np.all(np.abs(x - y) < 1e-14) and x.ann == y.ann

    def test_scw_invalid(self):
        a = acoustotreams.ScalarSphericalWaveBasis([[1, 0]])
        b = acoustotreams.ScalarCylindricalWaveBasis([[1, -1], [1, 1]])
        with pytest.raises(AttributeError):
            acoustotreams.rotate(1, basis=(a, b))

    def test_spw(self):
        b = acoustotreams.ScalarPlaneWaveBasisByUnitVector([[1, 0, 0], [0, 1, 0]])
        a = acoustotreams.ScalarPlaneWaveBasisByUnitVector(
            [[np.cos(2), np.sin(2), 0], [-np.sin(2), np.cos(2), 0]]
        )
        where = [True, False]
        x = acoustotreams.rotate(2, basis=b, where=where)
        y = acoustotreams.AcousticsArray([[1, 0], [0, 0]], basis=(a, b))
        assert np.all(np.abs(x - y) < 1e-14) and x.ann == y.ann

    def test_spwp(self):
        b = acoustotreams.ScalarPlaneWaveBasisByComp([[1, 0], [0, 1]])
        a = acoustotreams.ScalarPlaneWaveBasisByComp(
            [[np.cos(2), np.sin(2)], [-np.sin(2), np.cos(2)]]
        )
        where = [True, False]
        x = acoustotreams.rotate(2, basis=b, where=where)
        y = acoustotreams.AcousticsArray([[1, 0], [0, 0]], basis=(a, b))
        assert np.all(np.abs(x - y) < 1e-14) and x.ann == y.ann

    def test_matmul(self):
        assert (
            acoustotreams.Rotate(1, 2, 3)
            @ acoustotreams.AcousticsArray(np.eye(4), basis=acoustotreams.ScalarSphericalWaveBasis.default(1))
            == acoustotreams.rotate(1, 2, 3, basis=acoustotreams.ScalarSphericalWaveBasis.default(1))
        ).all()

    def test_rmatmul(self):
        assert (
            acoustotreams.AcousticsArray(np.eye(4), basis=acoustotreams.ScalarSphericalWaveBasis.default(1))
            @ acoustotreams.Rotate(1, 2, 3)
            == acoustotreams.rotate(1, 2, 3, basis=acoustotreams.ScalarSphericalWaveBasis.default(1))
        ).all()


class TestTranslate:
    def test_ssw(self):
        a = acoustotreams.ScalarSphericalWaveBasis([[1, 0]])
        b = acoustotreams.ScalarSphericalWaveBasis([[1, -1], [1, 1]])
        where = [True, False]
        x = acoustotreams.translate(
            [[0, 0, 0], [0, 1, 1]], k0=3, basis=(a, b), material=(1000, 1029, 0), where=where
        )
        y = acoustotreams.AcousticsArray(
            [
                [[acoustotreams.ssw.translate(1, 0, 1, -1, 0, 0, 0, singular=False), 0]],
                [
                    [
                        acoustotreams.ssw.translate(
                            1, 0, 1, -1, np.sqrt(2), np.pi / 4, np.pi / 2, singular=False
                        ),
                        0,
                    ]
                ],
            ],
            basis=(None, a, b),
            k0=(None, 3, 3),
            material=(None,) + 2 * (acoustotreams.AcousticMaterial(1000, 1029, 0),),
        )
        assert np.all(np.abs(x - y) < 1e-14) and x.ann == y.ann

    def test_ssw_invalid(self):
        a = acoustotreams.ScalarSphericalWaveBasis([[1, 0]])
        b = acoustotreams.ScalarCylindricalWaveBasis([[1, -1], [1, 1]])
        with pytest.raises(AttributeError):
            acoustotreams.translate([1, 0, 0], k0=1, basis=(b, a))

    def test_invalid_r(self):
        b = acoustotreams.ScalarSphericalWaveBasis([[1, -1], [1, 1]])
        with pytest.raises(ValueError):
            acoustotreams.translate([1, 0], basis=b, k0=1)

    def test_cw(self):
        a = acoustotreams.ScalarCylindricalWaveBasis([[0.1, 0]])
        b = acoustotreams.ScalarCylindricalWaveBasis([[0.1, -1], [0.1, 1]])
        where = [True, False]
        x = acoustotreams.translate(
            [[0, 0, 0], [0, 1, 1]], k0=3, basis=(a, b), material=(1000, 1029, 0), where=where
        )
        y = acoustotreams.AcousticsArray(
            [
                [
                    [
                        acoustotreams.scw.translate(
                            0.1, 0, 0.1, -1, 0, 0, 0, singular=False
                        ),
                        0,
                    ]
                ],
                [
                    [
                        acoustotreams.scw.translate(
                            0.1,
                            0,
                            0.1,
                            -1,
                            np.sqrt(1 - 0.01),
                            np.pi / 2,
                            1,
                            singular=False,
                        ),
                        0,
                    ]
                ],
            ],
            basis=(None, a, b),
            k0=(None, 3, 3),
            material=(None,) + 2 * (acoustotreams.AcousticMaterial(1000, 1029, 0),),
        )
        assert np.all(np.abs(x - y) < 1e-14) and x.ann == y.ann

    def test_spw(self):
        a = acoustotreams.ScalarPlaneWaveBasisByUnitVector([[1, 0, 0]])
        b = acoustotreams.ScalarPlaneWaveBasisByUnitVector([[1, 0, 0], [0.6, 0.8, 0]])
        where = [True, False]
        x = acoustotreams.translate([[0, 0, 0], [1, 1, 1]], k0=1, basis=(a, b), where=where)
        y = acoustotreams.AcousticsArray(
            [
                [[acoustotreams.spw.translate(1, 0, 0, 0, 0, 0), 0]],
                [[acoustotreams.spw.translate(1, 0, 0, 1, 1, 1), 0]],
            ],
            basis=(None, a, b),
            k0=(1, 1),
            material=((), ()),
        )
        assert np.all(np.abs(x - y) < 1e-14) and x.ann == y.ann

    def test_spwp(self):
        a = acoustotreams.ScalarPlaneWaveBasisByComp([[4, 0]])
        b = acoustotreams.ScalarPlaneWaveBasisByComp([[4, 0], [4, 1]])
        where = [True, False]
        x = acoustotreams.translate([[0, 0, 0], [0, 1, 1]], k0=5, basis=(a, b), where=where)
        y = acoustotreams.AcousticsArray(
            [
                [[acoustotreams.spw.translate(4, 0, 3, 0, 0, 0), 0]],
                [[acoustotreams.spw.translate(4, 0, 3, 0, 1, 1), 0]],
            ],
            basis=(None, a, b),
            k0=(None, 5, 5),
            material=(None, acoustotreams.AcousticMaterial(), acoustotreams.AcousticMaterial()),
            modetype=(None, "up", "up"),
        )
        assert np.all(np.abs(x - y) < 1e-14) and x.ann == y.ann


class TestExpand:
    def test_ssw_ssw_sing(self):
        a = acoustotreams.ScalarSphericalWaveBasis([[2, 0]], [0, 1, 1])
        b = acoustotreams.ScalarSphericalWaveBasis([[1, -1], [1, 1]])
        where = [True, False]
        x = acoustotreams.expand(
            (a, b), ("regular", "singular"), k0=3, material=(1000, 1029, 0), where=where
        )
        y = acoustotreams.AcousticsArray(
            [[acoustotreams.ssw.translate(2, 0, 1, -1, np.sqrt(2), 0.25 * np.pi, 0.5 * np.pi), 0]],
            basis=(a, b),
            k0=3,
            material=acoustotreams.AcousticMaterial(1000, 1029, 0),
            modetype=("regular", "singular"),
        )
        assert np.all(np.abs(x - y) < 1e-14) and x.ann == y.ann

    def test_scw_scw_sing(self):
        a = acoustotreams.ScalarCylindricalWaveBasis([[0.2, 0]], [0, 1, 1])
        b = acoustotreams.ScalarCylindricalWaveBasis([[0.2, -1], [0.2, 1]])
        where = [True, False]
        x = acoustotreams.expand(
            (a, b), ("regular", "singular"), k0=3, material=(1000, 1029, 0), where=where
        )
        y = acoustotreams.AcousticsArray(
            [
                [
                    acoustotreams.scw.translate(
                        0.2, 0, 0.2, -1, np.sqrt(1 - 0.04), np.pi / 2, 1
                    ),
                    0,
                ]
            ],
            basis=(a, b),
            k0=3,
            material=acoustotreams.AcousticMaterial(1000, 1029, 0),
            modetype=("regular", "singular"),
        )
        assert np.all(np.abs(x - y) < 1e-14) and x.ann == y.ann

    def test_ssw_scw(self):
        a = acoustotreams.ScalarSphericalWaveBasis([[1, 1]])
        b = acoustotreams.ScalarCylindricalWaveBasis([[0.3, 1], [0.1, 1]])
        where = [True, False]
        x = acoustotreams.expand((a, b), k0=3, material=(1000, 1029, 0), where=where)
        y = acoustotreams.AcousticsArray(
            [[acoustotreams.scw.to_ssw(1, 1, 0.3, 1, 1), 0]],
            basis=(a, b),
            k0=3,
            material=acoustotreams.AcousticMaterial(1000, 1029, 0),
            modetype="regular",
        )
        assert np.all(np.abs(x - y) < 1e-14) and x.ann == y.ann

    def test_spw_spw(self):
        a = acoustotreams.ScalarPlaneWaveBasisByUnitVector([[3, 0, 4], [0, 0, 5]])
        b = acoustotreams.ScalarPlaneWaveBasisByUnitVector([[3, 0, 4], [0, 5, 0]])
        where = [True, False]
        x = acoustotreams.expand((a, b), k0=2.5, material=(1000, 1200, 0), where=where)
        y = acoustotreams.AcousticsArray(
            [[1, 0], [0, 0]], basis=(a, b), k0=2.5, material=acoustotreams.AcousticMaterial(1000, 1200, 0),
        )
        assert np.all(np.abs(x - y) < 1e-14) and x.ann == y.ann

    def test_scw_spw(self):
        a = acoustotreams.ScalarCylindricalWaveBasis([[3, 1]], [1, 2, 3])
        b = acoustotreams.ScalarPlaneWaveBasisByUnitVector([[0, 4, 3], [0, 4, 3]])
        where = [True, False]
        x = acoustotreams.expand((a, b), k0=5, where=where)
        y = acoustotreams.AcousticsArray(
            [
                [
                    acoustotreams.spw.to_scw(3, 1, 0, 4, 3)
                    * acoustotreams.spw.translate(0, 4, 3, 1, 2, 3),
                    0,
                ]
            ],
            basis=(a, b),
            k0=5,
            material=acoustotreams.AcousticMaterial(),
            modetype=("regular", None),
        )
        assert np.all(np.abs(x - y) < 1e-14) and x.ann == y.ann

    def test_ssw_spw(self):
        a = acoustotreams.ScalarSphericalWaveBasis([[3, 1]], [1, 2, 3])
        b = acoustotreams.ScalarPlaneWaveBasisByUnitVector([[0, 4, 3]])
        where = [True, False]
        x = acoustotreams.expand((a, b), k0=5, where=where)
        y = acoustotreams.AcousticsArray(
            [
                [
                    acoustotreams.spw.to_ssw(3, 1, 0, 4, 3)
                    * acoustotreams.spw.translate(0, 4, 3, 1, 2, 3),
                    0,
                ]
            ],
            basis=(a, b),
            k0=5,
            material=acoustotreams.AcousticMaterial(),
            modetype=("regular", None),
        )
        assert np.all(np.abs(x - y) < 1e-14) and x.ann == y.ann


class TestExpandLattice:
    def test_ssw_1d(self):
        b = acoustotreams.ScalarSphericalWaveBasis([[1, -1], [1, 1]])
        where = [[True, False], [False, False]]
        lattice = acoustotreams.Lattice(1)
        x = acoustotreams.expandlattice(lattice, 0, basis=b, k0=3, where=where)
        y = acoustotreams.AcousticsArray(
            [
                [
                    acoustotreams.ssw.translate_periodic(
                        3, 0, lattice[...], [0, 0, 0], [[1], [-1]]
                    )[0, 0],
                    0,
                ],
                [0, 0],
            ],
            basis=b,
            k0=3,
            material=acoustotreams.AcousticMaterial(),
            modetype=("regular", "singular"),
            kpar=[np.nan, np.nan, 0],
            lattice=lattice,
        )
        assert np.all(np.abs(x - y) < 1e-14) and x.ann == y.ann

    def test_ssw_2d(self):
        b = acoustotreams.ScalarSphericalWaveBasis([[1, -1], [1, 1]])
        where = [[True, False], [False, False]]
        lattice = acoustotreams.Lattice([[1, 0], [0, 1]])
        x = acoustotreams.expandlattice(lattice, [0, 0], basis=b, k0=3, where=where)
        y = acoustotreams.AcousticsArray(
            [
                [
                    acoustotreams.ssw.translate_periodic(
                        3, [0, 0], lattice[...], [0, 0, 0], [[1], [-1]]
                    )[0, 0],
                    0,
                ],
                [0, 0],
            ],
            basis=b,
            k0=3,
            material=acoustotreams.AcousticMaterial(),
            modetype=("regular", "singular"),
            kpar=[0, 0, np.nan],
            lattice=lattice,
        )
        assert np.all(np.abs(x - y) < 1e-14) and x.ann == y.ann

    def test_ssw_3d(self):
        b = acoustotreams.ScalarSphericalWaveBasis([[1, -1], [1, 1]])
        where = [[True, False], [False, False]]
        lattice = acoustotreams.Lattice(np.eye(3))
        x = acoustotreams.expandlattice(lattice, [0, 0, 0], basis=b, k0=3, where=where)
        y = acoustotreams.AcousticsArray(
            [
                [
                    acoustotreams.ssw.translate_periodic(
                        3, [0, 0, 0], lattice[...], [0, 0, 0], [[1], [-1]]
                    )[0, 0],
                    0,
                ],
                [0, 0],
            ],
            basis=b,
            k0=3,
            material=acoustotreams.AcousticMaterial(),
            modetype=("regular", "singular"),
            kpar=[0, 0, 0],
            lattice=lattice,
        )
        assert np.all(np.abs(x - y) < 1e-14) and x.ann == y.ann

    def test_scw_ssw(self):
        a = acoustotreams.ScalarCylindricalWaveBasis([[0.3, 1]])
        b = acoustotreams.ScalarSphericalWaveBasis([[1, 1], [1, -1]])
        where = [True, False]
        x = acoustotreams.expandlattice(2, basis=(a, b), k0=3, material=(1000, 1029, 0), where=where)
        y = acoustotreams.AcousticsArray(
            [[acoustotreams.ssw.periodic_to_scw(0.3, 1, 1, 1, 1, 2), 0]],
            basis=(a, b),
            k0=3,
            material=acoustotreams.AcousticMaterial(1000, 1029, 0),
            modetype="singular",
            lattice=acoustotreams.Lattice(2),
            kpar=[np.nan, np.nan, 0.3],
        )
        assert np.all(np.abs(x - y) < 1e-14) and x.ann == y.ann

    def test_spw_ssw(self):
        a = acoustotreams.ScalarPlaneWaveBasisByUnitVector([[3, 0, 4]])
        b = acoustotreams.ScalarSphericalWaveBasis([[1, -1], [1, 1]])
        where = [True, False]
        lattice = acoustotreams.Lattice([[2, 0], [0, 2]])
        x = acoustotreams.expandlattice(
            lattice, [3, 0], basis=(a, b), k0=15, material=(1000, 1029, 0), where=where
        )
        y = acoustotreams.AcousticsArray(
            [[acoustotreams.ssw.periodic_to_spw(3, 0, 4, 1, -1, 4), 0]],
            basis=(a, b),
            k0=15,
            kpar=[3, 0, np.nan],
            lattice=lattice,
            material=acoustotreams.AcousticMaterial(1000, 1029, 0),
            modetype=(None, "singular"),
        )
        assert np.all(np.abs(x - y) < 1e-14) and x.ann == y.ann

    def test_scw_1d(self):
        b = acoustotreams.ScalarCylindricalWaveBasis([[0.1, -1], [0.1, 1]])
        where = [[True, False], [False, False]]
        lattice = acoustotreams.Lattice(1, "x")
        x = acoustotreams.expandlattice(lattice, 0, basis=b, k0=3, where=where)
        y = acoustotreams.AcousticsArray(
            [
                [
                    acoustotreams.scw.translate_periodic(
                        3, 0, lattice[...], [0, 0, 0], [[0], [0.1], [-1]]
                    )[0, 0],
                    0,
                ],
                [0, 0],
            ],
            basis=b,
            k0=3,
            material=acoustotreams.AcousticMaterial(),
            modetype=("regular", "singular"),
            kpar=[0, np.nan, 0.1],
            lattice=lattice,
        )
        assert np.all(np.abs(x - y) < 1e-14) and x.ann == y.ann

    def test_scw_2d(self):
        b = acoustotreams.ScalarCylindricalWaveBasis([[0.1, -1], [0.1, 1]])
        where = [[True, False], [False, False]]
        lattice = acoustotreams.Lattice([[1, 0], [0, 1]])
        x = acoustotreams.expandlattice(lattice, [0, 0], basis=b, k0=3, where=where)
        y = acoustotreams.AcousticsArray(
            [
                [
                    acoustotreams.scw.translate_periodic(
                        3, [0, 0], lattice[...], [0, 0, 0], [[0], [0.1], [-1]]
                    )[0, 0],
                    0,
                ],
                [0, 0],
            ],
            basis=b,
            k0=3,
            material=acoustotreams.AcousticMaterial(),
            modetype=("regular", "singular"),
            kpar=[0, 0, 0.1],
            lattice=lattice,
        )
        assert np.all(np.abs(x - y) < 1e-14) and x.ann == y.ann

    def test_spw_scw(self):
        a = acoustotreams.ScalarPlaneWaveBasisByUnitVector([[0, 3, 4]])
        b = acoustotreams.ScalarCylindricalWaveBasis([[4, -1], [4, 1]])
        where = [True, False]
        lattice = acoustotreams.Lattice(2, "x")
        x = acoustotreams.expandlattice(
            lattice, 3, basis=(a, b), k0=15, material=(1000, 1029, 0), where=where
        )
        y = acoustotreams.AcousticsArray(
            [[acoustotreams.scw.periodic_to_spw(0, 3, 4, 4, -1, 2), 0]],
            basis=(a, b),
            k0=15,
            material=acoustotreams.AcousticMaterial(1000, 1029, 0),
            modetype=(None, "singular"),
            lattice=lattice,
            kpar=[3, np.nan, np.nan],
        )
        assert np.all(np.abs(x - y) < 1e-14) and x.ann == y.ann


class TestPermute:
    def test_spw(self):
        a = acoustotreams.ScalarPlaneWaveBasisByUnitVector([[1, 2, 3]])
        b = acoustotreams.ScalarPlaneWaveBasisByUnitVector([[2, 3, 1]])
        assert acoustotreams.permute(basis=b).basis[0] == a

    def test_spwp(self):
        a = acoustotreams.ScalarPlaneWaveBasisByComp([[1, 2]], "yz")
        b = acoustotreams.ScalarPlaneWaveBasisByComp([[1, 2]])
        assert acoustotreams.permute(basis=b, k0=5, material=1).basis[0] == a


class TestPField:
    def test_ssw_r(self):
        modes = [[0, 3, -2], [1, 1, 1]]
        positions = np.array([[0, 0, 0], [1, 0, 0]])
        b = acoustotreams.ScalarSphericalWaveBasis(modes, positions)
        r = np.array([[0, 1, 2], [3, 4, 5]])
        k0 = 4
        material = (1000, 686, 0)
        x = acoustotreams.pfield(r, basis=b, k0=k0, material=material)
        rsph = acoustotreams.car2sph(r[:, None] - positions)
        y = acoustotreams.ssw_rPsi(
            [3, 1],
            [-2, 1],
            k0 * rsph[..., 0] * 0.5,
            rsph[..., 1],
            rsph[..., 2],
        )
        y = np.array([y]).T
        assert np.all(np.abs(y.swapaxes(-1, -2) - x) < 1e-14)

    def test_ssw_s(self):
        modes = [[0, 3, -2], [1, 1, 1]]
        positions = np.array([[0, 0, 0], [1, 0, 0]])
        b = acoustotreams.ScalarSphericalWaveBasis(modes, positions)
        r = np.array([[0, 1, 2], [3, 4, 5]])
        k0 = 4
        material = (1000, 686, 0)
        x = acoustotreams.pfield(
            r,
            basis=b,
            k0=k0,
            material=material,
            modetype="singular",
        )
        rsph = acoustotreams.car2sph(r[:, None] - positions)
        y = acoustotreams.ssw_Psi(
            [3, 1],
            [-2, 1],
            k0 * rsph[..., 0] * 0.5,
            rsph[..., 1],
            rsph[..., 2],
        )
        y = np.array([y]).T
        assert np.all(np.abs(y.swapaxes(-1, -2) - x) < 1e-14)

    def test_scw_r(self):
        modes = [[0, 0.3, -2], [1, 0.1, 1]]
        positions = np.array([[0, 0, 0], [1, 0, 0]])
        b = acoustotreams.ScalarCylindricalWaveBasis(modes, positions)
        r = np.array([[0, 1, 2], [3, 4, 5]])
        k0 = 4
        material = (1000, 686, 0)
        x = acoustotreams.pfield(r, basis=b, k0=k0, material=material)
        rcyl = acoustotreams.car2cyl(r[:, None] - positions)
        y = acoustotreams.scw_rPsi(
            [0.3, 0.1],
            [-2, 1],
            rcyl[..., 0] * [np.sqrt(4 - 0.3 ** 2), np.sqrt(4 - 0.1 ** 2)],
            rcyl[..., 1],
            rcyl[..., 2],
        )
        y = np.array([y]).T
        assert np.all(np.abs(y.swapaxes(-1, -2) - x) < 1e-14)

    def test_scw_s(self):
        modes = [[0, 0.3, -2], [1, 0.1, 1]]
        positions = np.array([[0, 0, 0], [1, 0, 0]])
        b = acoustotreams.ScalarCylindricalWaveBasis(modes, positions)
        r = np.array([[0, 1, 2], [3, 4, 5]])
        k0 = 4
        material = (1000, 686, 0)
        x = acoustotreams.pfield(
            r,
            basis=b,
            k0=k0,
            material=material,
            modetype="singular",
        )
        rcyl = acoustotreams.car2cyl(r[:, None] - positions)
        y = acoustotreams.scw_Psi(
            [0.3, 0.1],
            [-2, 1],
            rcyl[..., 0] * [np.sqrt(4 - 0.3 ** 2), np.sqrt(4 - 0.1 ** 2)],
            rcyl[..., 1],
            rcyl[..., 2],
        )
        y = np.array([y]).T
        assert np.all(np.abs(y.swapaxes(-1, -2) - x) < 1e-14)

    def test_spw(self):
        modes = [[0, 3, 4], [0, 4, 3]]
        b = acoustotreams.ScalarPlaneWaveBasisByUnitVector(modes)
        r = np.array([[0, 1, 2], [3, 4, 5]])
        x = acoustotreams.pfield(r, basis=b, k0=5)
        y = acoustotreams.spw_Psi(
            [0, 0],
            [3, 4],
            [4, 3],
            r[..., None, 0],
            r[..., None, 1],
            r[..., None, 2],
        )
        y = np.array([y]).T
        assert np.all(np.abs(y.swapaxes(-1, -2) - x) < 1e-14)

    def test_spwp(self):
        modes = [[0, 3], [0, 4]]
        b = acoustotreams.ScalarPlaneWaveBasisByComp(modes)
        r = np.array([[0, 1, 2], [3, 4, 5]])
        x = acoustotreams.pfield(r, basis=b, k0=5)
        y = acoustotreams.spw_Psi(
            [0, 0],
            [3, 4],
            [4, 3],
            r[..., None, 0],
            r[..., None, 1],
            r[..., None, 2],
        )
        y = np.array([y]).T
        assert np.all(np.abs(y.swapaxes(-1, -2) - x) < 1e-14)