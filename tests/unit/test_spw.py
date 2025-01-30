import numpy as np

from acoustotreams import spw


def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


class TestToScw:
    def test(self):
        assert isclose(spw.to_scw(1, 4, 3, 2, 1), np.power((2 + 3j), 4) / 169)

    def test_complex(self):
        assert isclose(spw.to_scw(1, 4, 3, 2j, 1), 25)

    def test_one(self):
        assert spw.to_scw(1, 0, 3, 2, 1) == 1


class TestToSsw:
    def test(self):
        assert isclose(
            spw.to_ssw(5, 4, 3, 2, 1), 3.0179547476477375 - 2.9928051247506713j
        )


class TestTranslate:
    def test(self):
        assert spw.translate(1, 2, 3, 4, 5, 6) == np.exp(32j)


class TestPermuteXyz:
    def test(self):
        assert isclose(
            spw.permute_xyz(1, 2, 3), 1
        )

    def test_inv(self):
        assert isclose(
            spw.permute_xyz(1, 2, 1j, inverse=True), 1
        )