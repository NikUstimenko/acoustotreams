import numpy as np

from acoustotreams import ssw

def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


class TestPeriodicToScw:
    def test(self):
        assert isclose(
            ssw.periodic_to_scw(1, 2, 3, 2, 3, 2),
            -0.15855121307843276j,
        )

    def test_zero(self):
        assert ssw.periodic_to_scw(1, 2, 3, -2, 3, 2) == 0j


class TestPeriodicToSpw:
    def test(self):
        assert isclose(
            ssw.periodic_to_spw(1, 2, 3 + 4j, 5, -4, 2),
            -0.000283588647059046 + 0.008247703151966499j,
        )

    def test_zero(self):
        assert ssw.periodic_to_spw(1, 2, 0, 5, -4, 2) == 0j


    def test_kz_neg(self):
        assert isclose(
            ssw.periodic_to_spw(1, 2, -3, 5, -4, 2),
            -0.04032911669459158 - 0.011762659035922562j,
        )


class TestRotate:
    def test(self):
        assert isclose(
            ssw.rotate(6, 5, 6, 4, 3, 2, 1), 0.04974201284017233 - 0.0075403653936369315j
        )
    def test_zero(self):
        assert ssw.rotate(8, 7, 6, 5, 4, 3, 2) == 0j


class TestTranslate:
    def test_s_real(self):
        assert isclose(
            ssw.translate(5, 4, 3, 2, 8, 7, 6),
            -0.021827556564340267 - 0.19541673271082188j,
        )

    def test_r_real(self):
        assert isclose(
            ssw.translate(5, 4, 3, 2, 8, 7, 6, singular=False),
            -0.10402575627867332 - 0.0661458099663242j,
        )

    def test_s_kr_zero(self):
        assert isclose(ssw.translate(5, 4, 3, 2, 1e-30j, 7, 6), 0j)

    def test_r_kr_zero(self):
        assert isclose(ssw.translate(5, 4, 3, 2, 1e-30j, 7, 6, singular=False), 0j)


class TestTranslatePeriodic:
    def test_0(self):
        assert isclose(
            ssw.translate_periodic(
                1, 0, 1, [0, 0, 0], ([1], [0]), in_=([1], [0])
            )[0, 0],
            -1. + 17.29826864j
        )

    def test_1(self):
        assert isclose(
            ssw.translate_periodic(
                1, [0, 0], [[1, 0], [0, 1]], [0, 0, 0], ([2], [-1])
            )[0, 0],
            -1. + 381.01677719j,
        )

    def test_2(self):
        assert isclose(
            ssw.translate_periodic(
                1,
                [0, 0, 0],
                [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
                [0, 0, 0],
                ([3], [3])[0, 0]
            ),
            -1. - 1203.51273971j,
        )

    def test_3(self):
        assert (
            ssw.translate_periodic(
                1,
                [0, 0],
                [[1, 0], [0, 1]],
                [0, 0, 0],
                ([1], [0]),
                in_=([2], [0]),
            )[0, 0]
            == 0. + 0j
        )
