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

    def test_coreshell_ff(self):
        expect = np.array(
            [
                -0.008088738845585024 - 0.0895729375954197j,
                -1.0441230501976207e-06 - 0.0010218228613632786j,
                -5.479473213464122e-08 + 0.0002340827399278665j
            ]
        )
        assert np.all(
            np.abs(
                cf.mie_acoustics([0, 1, 2], [1, 2], [800, 1000, 900], [1000, 1200, 800], [0 , 0, 0])
                - expect
            )
            < EPSSQ
        )

    def test_coreshell_sf(self):
        expect = np.array(
            [
                -0.006854218938361541 - 0.08250599142551134j,
                -8.62868748368375e-07 - 0.000928906886522314j,
                -2.350578936408102e-07 + 0.0004848276378141699j
            ]
        )
        assert np.all(
            np.abs(
                cf.mie_acoustics([0, 1, 2], [1, 2], [800, 1000, 900], [1000, 1200, 800], [500, 0, 0])
                - expect
            )
            < EPSSQ
        )

    def test_coreshell_ss(self):
        expect = np.array(
            [
                -0.0034154783440397945 - 0.05834220472112098j,
                -1.514571694677478e-06 + 0.0012306784310913421j,
                -9.607196474592923e-05 + 0.00980115987643904j
            ]
        )
        assert np.all(
            np.abs(
                cf.mie_acoustics([0, 1, 2], [1, 2], [800, 1000, 900], [1000, 1200, 800], [500, 600, 0])
                - expect
            )
            < EPSSQ
        )

    def test_coreshell_fs(self):
        expect = np.array(
            [
                -0.004115607551091826 - 0.06402085070957136j,
                -1.3567878242366934e-06 + 0.001164811565602216j,
                -0.0002426314464083701 + 0.015574741615501845j
            ]
        )
        assert np.all(
            np.abs(
                cf.mie_acoustics([0, 1, 2], [1, 2], [800, 1000, 900], [1000, 1200, 800], [0, 600, 0])
                - expect
            )
            < EPSSQ
        )

    def test_coreshell_softf(self):
        expect = np.array(
            [
                -0.5718522069566513 - 0.4948103276564248j,
                -0.0219397833519953 - 0.14648695934608927j,
                -7.317407093548152e-05 - 0.008553871432914141j
            ]
        )
        assert np.all(
            np.abs(
                cf.mie_acoustics([0, 1, 2], 2, [0, 900], [0, 800], [0, 0])
                - expect
            )
            < EPSSQ
        )

    def test_coreshell_hardf(self):
        expect = np.array(
            [
                -0.0219397833519953 - 0.14648695934608927j,
                -0.006054311374897304 + 0.07757355663286992j,
                -3.545724003859597e-05 + 0.005954492658717705j
            ]
        )
        assert np.all(
            np.abs(
                cf.mie_acoustics([0, 1, 2], 2, [np.inf, 900], [0, 800], [0, 0])
                - expect
            )
            < EPSSQ
        )



class TestMieCyl:
    def test_real(self):
        expect = 1.7330362995685293e-20 + 0.0002685792010794467j
        assert (
            np.abs(
                cf.mie_acoustics_cyl(0.5, -2, 1, 2, [1000, 900], [1200, 800], [0, 0])
                - expect
            )
            < EPSSQ
        )

    def test_complex(self):
        expect = -0.00038763204996962404 + 0.0002836447825423364j
        assert (np.abs(
                    cf.mie_acoustics_cyl(
                        0.5, 
                        -2,
                        1,
                        2, 
                        [1000 + 100j, 900], 
                        [1200 - 150j, 800], 
                        [0, 0]
                        )
                    - expect
                ) 
                < EPSSQ
        )


class TestFresnel:
    def test_real(self):
        expect = [
            [
                [[0.37522557959686687]],
                [[0.6247744204031331]],
            ],
            [
                [[-0.6247744204031331]],
                [[1.6247744204031331]],
            ],
        ]
        assert np.all(
            np.abs(
                cf.fresnel_acoustics(
                    [np.sqrt(8), np.sqrt(24)],
                    [1000, 400],
                )
                - expect
            )
            < EPSSQ
        )

    def test_evanescent(self):
        expect = [
            [
                [[0.10126582278481014 + 0.43849387533389306j]],
                [[0.89873417721519 - 0.43849387533389306j]],
            ],
            [
                [[-0.89873417721519 + 0.43849387533389306j]],
                [[1.89873417721519 - 0.43849387533389306j]],
            ],
        ]
        assert np.all(
            np.abs(
                cf.fresnel_acoustics(
                    [1j * np.sqrt(8), np.sqrt(24)],
                    [1000, 400],
                )
                - expect
            )
            < EPSSQ
        )


    def test_complex(self):
        expect = [
            [
                [[0.3960187651651478 + 0.07003355217085547j]],
                [[0.6039812348348521 - 0.07003355217085547j]],
            ],
            [
                [[-0.6039812348348521 + 0.07003355217085547j]],
                [[1.603981234834852 - 0.07003355217085547j]],
            ],
        ]
        assert np.all(
            np.abs(
                cf.fresnel_acoustics(
                    [np.sqrt(8), np.sqrt(24)],
                    [1000 + 250j, 400 + 200j],
                )
                - expect
            )
            < EPSSQ
        )


