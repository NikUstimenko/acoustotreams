import numpy as np
import scipy.special as sc

from acoustotreams import scw


def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)

class TestPeriodicToSpw:
    def test(self):
        assert isclose(
            scw.periodic_to_spw(6, -5, 4, 4, 3, 2),
            2 * np.power((-5 - 6j) / np.sqrt(61), 3) / 10,
        )

    def test_ky_zero(self):
        assert isclose(scw.periodic_to_spw(6, 0, 4, 4, 3, 2), 5e19 + 5e19j)


class TestRotate:
    def test(self):
        assert isclose(scw.rotate(3, 2, 3, 2, 4), np.exp(-8j))

    def test_zero(self):
        assert scw.rotate(3, 2, 2, 2, 4) == 0


class TestToSsw:
    def test(self):
        assert isclose(scw.to_ssw(4, 3, 2, 3, 5), -4.843721204315066j)

    def test_zero(self):
        assert scw.to_ssw(4, 3, 2, 2, 5) == 0j


class TestTranslate:
    def test_s(self):
        assert isclose(
            scw.translate(3, 2, 3, -2, 4, 5, 6), sc.hankel1(-4, 4) * np.exp(-2j)
        )

    def test_s_opposite(self):
        assert scw.translate(3, 2, 0, 2, 4 + 1j, 5, 6) == 0j

    def test_s_zero(self):
        assert scw.translate(3, 2, 0, 2, 0, 5, 0) == 0j

    def test_r(self):
        assert isclose(
            scw.translate(3, 2, 3, -2, 4, 5, 6, singular=False),
            sc.jv(-4, 4) * np.exp(-2j),
        )

    def test_r_opposite(self):
        assert scw.translate(3, 2, 0, 2, 4 + 1j, 5, 6, singular=False) == 0j


class TestTranslatePeriodic:
    def test(self):
        assert isclose(
            scw.translate_periodic(1, 0, 2, [0, 0, 0], ([0], [0]))[0, 0],
            0.7610235793674727j,
        )

    def test_ks(self):
        assert isclose(
            scw.translate_periodic(2, 0, 2, [0, 0, 0], ([1], [3]))[0, 0],
            -0.42264973081037416 + 0.30599871466756257j,
        )

    def test_ks_same(self):
        assert isclose(
            scw.translate_periodic(
                1, [0, 1], [[2, 0], [0, 2]], [0, 0, 0], ([0], [0])
            )[0, 0],
            -1. + 0.28364898059265115j,
        )