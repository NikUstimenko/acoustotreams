import copy

import numpy as np
import pytest

import acoustotreams
from acoustotreams import AcousticTMatrix


def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


class TestInit:
    def test_simple(self):
        tm = AcousticTMatrix(np.eye(4), k0=1)
        assert (
            np.all(tm == np.eye(4))
            and tm.k0 == 1
            and tm.material == (1.3, 343., 0.)
            and tm.basis == acoustotreams.ScalarSphericalWaveBasis.default(1)
        )

    def test_complex(self):
        tm = AcousticTMatrix(
            np.diag([1, 2]),
            k0=3,
            material=[200, 800, 100],
            basis=acoustotreams.ScalarSphericalWaveBasis([[0, 0], [1, 0]], [1, 0, 0]),
        )
        assert (
            np.all(tm == np.diag([1, 2]))
            and tm.k0 == 3
            and tm.material == (200., 800., 100.)
            and tm.basis == acoustotreams.ScalarSphericalWaveBasis([[0, 0], [1, 0]], [1, 0, 0])
        )

class TestSphere:
    def test(self):
        tm = AcousticTMatrix.sphere(2, 3, 4, [(200, 1000, 500), (900, 800, 0)])
        del tm.modetype
        m = acoustotreams.coeffs.mie_acoustics([0, 1, 2], [12], [200, 900], [1000, 800], [500, 0])
        assert (
            tm[0, 0] == m[0]
            and np.all(np.diag(tm)[1:4] == m[1])
            and np.all(np.diag(tm)[4::] == m[2])
            and np.all(np.diag(tm, -1) == 0.+0.j)
            and np.all(np.diag(tm, 1) == 0.+0.j)
            and np.all(np.triu(tm, 2) == 0.+0.j)
            and np.all(np.tril(tm, -2) == 0.+0.j)
            and tm.k0 == 3
            and tm.material == (900, 800, 0)
            and tm.basis == acoustotreams.ScalarSphericalWaveBasis.default(2)
        )

    def test_invalid_material(self):
        with pytest.raises(ValueError):
            AcousticTMatrix.sphere(2, 3, 4, [(200, 1000, 500), (900, 800, 10)])


class TestProperties:
    def test_xs_ext_avg(self):
        tm = AcousticTMatrix.sphere(2, 3, [4], [(200 + 10j, 1000 - 100j, 500 - 50j), (900, 800, 0)])
        assert isclose(tm.xs_ext_avg, 24.60384476871747)

    def test_xs_ext_avg_ct_zero(self):
        tm = AcousticTMatrix.sphere(2, 3, [4], [(200 + 10j, 1000 - 100j, 0), (900, 800, 0)])
        assert isclose(tm.xs_ext_avg, 25.953193306013617)

    def test_xs_sca_avg(self):
        tm = AcousticTMatrix.sphere(2, 3, [4], [(200 + 10j, 1000 - 100j, 500 - 50j), (900, 800, 0)])
        assert isclose(tm.xs_sca_avg, 15.502943925268221)

    def test_xs_sca_avg_ct_zero(self):
        tm = AcousticTMatrix.sphere(2, 3, [4], [(200 + 10j, 1000 - 100j, 0), (900, 800, 0)])
        assert isclose(tm.xs_sca_avg, 15.258304851607877)

    def test_modes(self):
        tm = AcousticTMatrix.sphere(2, 3, [4], [(200 + 10j, 1000 - 100j, 500 - 50j), (900, 800, 0)])
        l, m = tm.basis.lm
        assert (
            np.all(l == [0] + 3 * [1] + 5 * [2])
            and np.all(m == [0, -1, 0, 1, -2, -1, 0, 1, 2])
        )


class TestXs:
    def test(self):
        tm = AcousticTMatrix.sphere(2, 3, [4], [(200 + 10j, 1000 - 100j, 500 - 50j), (900, 800, 0)])
        inc = acoustotreams.plane_wave_scalar([0, 0, 1], k0=tm.k0, material=tm.material)
        xs = tm.xs(inc, 0.125)
        assert isclose(xs[0], 62.0117757010729,) and isclose(xs[1], 98.41537907486989)


class TestSca:
    def test(self):
        tm = AcousticTMatrix.sphere(1, 3, [4], [(200 + 10j, 1000 - 100j, 500 - 50j), (900, 800, 0)])
        inc = acoustotreams.plane_wave_scalar([0, 0, 1], k0=tm.k0, material=tm.material)
        sca = tm.sca(inc)
        assert (isclose(sca[0], -2.7586039218306673 + 0.17803571960546696j,) 
                and isclose(sca[2], 2.4084972606519726 - 2.410569039823313j)
                and sca[1] == sca[3] == 0. + 0.j)


class TestTranslate:
    def test(self):
        tm = AcousticTMatrix.sphere(3, 0.1, [0.2], [(200 + 10j, 1000 - 100j, 500 - 50j), (900, 800, 0)])
        m = copy.deepcopy(tm)
        rs = np.array([[0.1, 0.2, 0.3], [-0.4, -0.5, -0.4]])
        tm = tm.translate(rs[0])
        tm = tm.translate(rs[1])
        tm = tm.translate(-rs[0] - rs[1])
        assert np.all(np.abs(tm - m) < 1e-8)

    def test_ct_zero(self):
        tm = AcousticTMatrix.sphere(3, 0.1, [0.2], [(200 + 10j, 1000 - 100j, 0), (900, 800, 0)])
        m = copy.deepcopy(tm)
        rs = np.array([[0.1, 0.2, 0.3], [-0.4, -0.5, -0.4]])
        tm = tm.translate(rs[0])
        tm = tm.translate(rs[1])
        tm = tm.translate(-rs[0] - rs[1])
        assert np.all(np.abs(tm - m) < 1e-8)


class TestClusterRotate:
    def test(self):
        tms = [AcousticTMatrix.sphere(3, 0.1, [0.1], [i * i * 10, 1000]) for i in range(1, 5)]
        rs1 = np.array([[0, 0, 0], [0.2, 0, 0], [0, 0.2, 0], [0, 0, 0.2]])
        tm1 = AcousticTMatrix.cluster(tms, rs1)
        tm1 = tm1.interaction.solve().expand(acoustotreams.ScalarSphericalWaveBasis.default(3))
        tm1 = tm1.rotate(1, 2, 3)
        a = np.array([[np.cos(1), -np.sin(1), 0], [np.sin(1), np.cos(1), 0], [0, 0, 1]])
        b = np.array([[np.cos(2), 0, np.sin(2)], [0, 1, 0], [-np.sin(2), 0, np.cos(2)]])
        c = np.array([[np.cos(3), -np.sin(3), 0], [np.sin(3), np.cos(3), 0], [0, 0, 1]])
        rs2 = (a @ b @ c @ rs1.T).T
        tm2 = AcousticTMatrix.cluster(tms, rs2)
        tm2 = tm2.interaction.solve().expand(acoustotreams.ScalarSphericalWaveBasis.default(3))
        assert np.all(np.abs(tm1 - tm2) < 1e-16)