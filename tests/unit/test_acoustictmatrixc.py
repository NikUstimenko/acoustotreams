import copy

import numpy as np
import pytest

import acoustotreams
from acoustotreams import AcousticTMatrixC


def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


class TestInit:
    def test_simple(self):
        tm = AcousticTMatrixC(np.eye(3), k0=1)
        assert (
            np.all(tm == np.eye(3))
            and tm.k0 == 1
            and tm.material.rho == 1.3
            and tm.material.c == 343.
            and tm.material.ct == 0.
            and np.all(tm.basis.positions == [[0, 0, 0]])
            and np.all(tm.basis.kz == 3 * [0])
            and np.all(tm.basis.m == [-1, 0, 1])
            and np.all(tm.basis.pidx == [0, 0, 0])
            and tm.ks == 1
        )

    def test_complex(self):
        tm = AcousticTMatrixC(
            np.diag([1, 2]),
            k0=3,
            material=acoustotreams.AcousticMaterial(900, 686, 0),
            basis=acoustotreams.ScalarCylindricalWaveBasis([[1, -1], [1, 0]], [1, 0, 0]),
        )
        assert (
            np.all(tm == np.diag([1, 2]))
            and tm.k0 == 3
            and tm.material.rho == 900
            and tm.material.c == 686
            and tm.material.ct == 0
            and np.all(tm.basis.positions == [[1, 0, 0]])
            and np.all(tm.basis.kz == [1, 1])
            and np.all(tm.basis.m == [-1, 0])
            and np.all(tm.basis.pidx == [0, 0])
            and tm.ks == 1.5
        )


class TestCylinder:
    def test(self):
        tm = AcousticTMatrixC.cylinder([1], 2, 3, 4, [(200, 1000, 500), (900, 686, 0)])
        m = acoustotreams.coeffs.mie_cyl_acoustics(1, [-2, -1, 0, 1, 2], 3, [4], [200, 900], [1000, 686], [500, 0])
        assert (
            tm[0, 0] == m[0]
            and tm[1, 1] == m[1]
            and tm[2, 2] == m[2]
            and tm[3, 3] == m[3]
            and tm[4, 4] == m[4]
            and tm.k0 == 3
            and tm.material.rho == 900
            and tm.material.c == 686
            and tm.material.ct == 0
            and np.all(tm.basis.positions == [[0, 0, 0]])
            and np.all(tm.basis.kz == 5 * [1])
            and np.all(tm.basis.m == [-2, -1, 0, 1, 2])
            and np.all(tm.basis.pidx == 5 * [0])
            and tm.ks == 1.5
        )


class TestProperties:
    def test_xw_ext_avg(self):
        tm = AcousticTMatrixC.cylinder([-1, 1], 1, 3, [4], [(200 + 10j, 1000 - 100j, 500), (900, 686, 0)])
        assert isclose(tm.xw_ext_avg, 2.8575610823570337)

    def test_xw_ext_avg_ct_zero(self):
        tm = AcousticTMatrixC.cylinder([-1, 1], 1, 3, [4], [(200 + 10j, 1000 - 100j, 0), (900, 686, 0)])
        assert isclose(tm.xw_ext_avg, 2.976939980231719)

    def test_xw_sca_avg(self):
        tm = AcousticTMatrixC.cylinder([-1, 1], 1, 3, [4], [(200 + 10j, 1000 - 100j, 500), (900, 686, 0)])
        assert isclose(tm.xw_sca_avg, 2.5126294545531063)

    def test_xw_sca_avg_ct_zero(self):
        tm = AcousticTMatrixC.cylinder([-1, 1], 1, 3, [4], [(200 + 10j, 1000 - 100j, 0), (900, 686, 0)])
        assert isclose(tm.xw_sca_avg, 1.3598694262205515)

    def test_krho(self):
        tm = AcousticTMatrixC.cylinder([0, 5], 1, 3, [1], [(200 + 10j, 1000 - 100j,), ()])
        assert np.all(
            np.abs(tm.krhos - [3, 3, 3, 4j, 4j, 4j]) < 1e-16
        )

    def test_modes(self):
        tm =  AcousticTMatrixC.cylinder([-1, 1], 1, 3, [4], [(200 + 10j, 1000 - 100j,), (900, 686,)])
        kz, m = tm.basis.zm
        assert (
            np.all(kz == 3 * [-1] + 3 * [1])
            and np.all(m == [-1, 0, 1, -1, 0, 1])
        )

    def test_fullmodes(self):
        tm = AcousticTMatrixC.cylinder([-1, 1], 1, 3, [4], [(200 + 10j, 1000 - 100j,), (900, 686,)])
        pidx, kz, m = tm.basis[()]
        assert (
            np.all(kz == 3 * [-1] + 3 * [1])
            and np.all(m == [-1, 0, 1, -1, 0, 1])
            and np.all(pidx == 6 * [0])
        )


class TestXw:
    def test(self):
        kz = 1
        tm = AcousticTMatrixC.cylinder(kz, 1, 3, [4], [(200 + 10j, 1000 - 100j,), ()])
        inc = acoustotreams.plane_wave_scalar([np.sqrt(tm.k0 * tm.k0 - kz * kz), 0, kz], k0=tm.k0, material=tm.material)
        xw = tm.xw(inc, 0.1255)
        assert isclose(xw[0], 5.892670869743773,) and isclose(xw[1], 5.902641304396861)


class TestSca:
    def test(self):
        kz = 1
        tm = AcousticTMatrixC.cylinder(kz, 1, 3, [4], [(200 + 10j, 1000 - 100j,), ()])
        inc = acoustotreams.plane_wave_scalar([np.sqrt(tm.k0 * tm.k0 - kz * kz), 0, kz], k0=tm.k0, material=tm.material)
        sca = tm.sca(inc)
        assert (isclose(sca[0], 0.3506619987521309 + 0.14429243052400523j,) 
                and isclose(sca[1], -0.8225873645046986 - 0.38088979145646645j)
                and isclose(sca[2], -0.3506619987521309 - 0.14429243052400523j))