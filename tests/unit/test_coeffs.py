import numpy as np

import acoustotreams.coeffs as cf

EPS = 2e-7
EPSSQ = 4e-14


class TestMie:
    def test_real(self):
        expect = np.array(
            [
                -0.00616101729970458 - 0.07824997869352635j,
                -2.968300985715852e-06 + 0.0017228732323954813j,
                -0.00011642422284159261 + 0.010789377565083367j
            ]
        )
        assert np.all(
            np.abs(
                cf.mie_acoustics([0, 1, 2], 2, [1000, 900], [1200, 800], [500, 0])
                - expect
            )
            < EPSSQ
        )

    def test_complex(self):
        expect = np.array(
            [
                -0.018834091697937875 - 0.07924217855202398j,
                -0.00629117938005657 + 0.0018032888761669984j,
                -0.0011494426526906401 + 0.010442610865723645j
            ]
        )
        assert np.all(
            np.abs(
                cf.mie_acoustics(
                    [0, 1, 2], 
                    2, 
                    [1000 + 100j, 900], 
                    [1200 - 150j, 800], 
                    [500 - 50j, 0]
                    )
                - expect
            )
            < EPSSQ
        )


class TestMieCyl:
    def test_real(self):
        expect = -0.0003022791824316385 + 0.011902814173445683j
        assert (
            np.abs(
                cf.mie_acoustics_cyl(0.5, -2, 1, 2, [1000, 900], [1200, 800], [500, 0])
                - expect
            )
            < EPSSQ
        )

    def test_complex(self):
        expect = 0.00037060887669232586 + 0.012983414050305396j
        assert (np.abs(
                    cf.mie_acoustics_cyl(
                        0.5, 
                        -2,
                        1,
                        2, 
                        [1000 + 100j, 900], 
                        [1200 - 150j, 800], 
                        [500 - 50j, 0]
                        )
                    - expect
                ) 
                < EPSSQ
        )
