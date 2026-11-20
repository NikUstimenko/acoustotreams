import numpy as np
import pytest
import acoustotreams
from acoustotreams import AcousticSMatrices


def test_init():
    b = acoustotreams.ScalarPlaneWaveBasisByComp.default([0, 0])
    sm = AcousticSMatrices(np.zeros((2, 2, 1, 1)), basis=b, k0=1)
    assert (sm["up", "down"] == np.zeros((1, 1))).all()


def test_interface():
    b = acoustotreams.ScalarPlaneWaveBasisByComp.default([1, 2])
    sm = AcousticSMatrices.interface(b, 12, [(200, 1029,), (900, 686,)])
    m = acoustotreams.coeffs.fresnel_acoustics(
        [np.sqrt(11), np.sqrt(31)], [200, 900]
    )
    assert (
        sm[0, 0] == m[0, 0]
        and sm[0, 1] == m[0, 1]
        and sm[1, 0] == m[1, 0]
        and sm[1, 1] == m[1, 1]
        and sm.k0 == 12
        and sm.basis == b
        and sm.material == ((900, 686,), (200, 1029,))
    )


def test_slab():
    b = acoustotreams.ScalarPlaneWaveBasisByComp.default([1, 2])
    sm = AcousticSMatrices.slab(3, b, 12, [(200, 1029,), (900, 686,), (100, 343,)])
    stest = AcousticSMatrices.interface(b, 12, [(200, 1029,), (900, 686,)])
    stest = stest.add(AcousticSMatrices.propagation([0, 0, 3], b, 12, (900, 686,)))
    stest = stest.add(AcousticSMatrices.interface(b, 12, [(900, 686,), (100, 343,)]))
    assert sm == stest


class TestArray:
    basis = acoustotreams.ScalarPlaneWaveBasisByComp.diffr_orders([0.5, -1], np.eye(2), 7)
    expect = np.zeros((2, 2, 5, 5), complex)
    expect[0, 0, :, :] = [
        [
            0.979645694563844 + 0.03743860913339489j,
            -0.017857473673999053 - 0.005077635955552676j,
            -0.013734573631049728 - 0.0046245525883680996j,
            -0.024443469729725567 - 0.028669459846426687j,
            -0.018376164691808704 - 0.0016794296265635839j,
        ],
        [
            -0.019823605196670754 - 0.0056366904049734885j,
            0.9708849843014916 + 0.04031767197824951j,
            -0.0069598869114165084 - 0.05888988234327362j,
            -0.028356227993153422 - 0.05907316358840312j,
            -0.018795029420162816 - 0.031674557647168j,
        ],
        [
            -0.017230073084042043 - 0.005801518213746545j,
            -0.00786523130953012 - 0.06655029777295496j,
            0.9733774204079847 + 0.04866102775883646j,
            -0.024748282499230397 - 0.0621489255327363j,
            -0.016937172015851695 - 0.032042486369294125j,
        ],
        [
            -0.029654528360355002 - 0.034781449585204335j,
            -0.030989459802982112 - 0.06455884854994856j,
            -0.02393322891270318 - 0.06010212876387582j,
            0.9374235611255546 + 0.004370404478317288j,
            -0.021882022410520837 - 0.09240662320571531j,
        ],
        [
            -0.02097622301763208 - 0.0019170534755230829j,
            -0.019326481970050186 - 0.0325701946824621j,
            -0.015411374639662323 - 0.029155915837738133j,
            -0.020588831922819812 - 0.08694554817853925j,
            0.9720874797759053 + 0.05307586568278426j,
        ],
    ]
    expect[0, 1, :, :] = [
        [
            -0.011258679597676765 + 0.022395921742739886j,
            -0.014099858804436134 + 0.020625426955153875j,
            -0.00718927006735804 + 0.010009695657248845j,
            -0.01777032615770584 - 0.011935898167403364j,
            -0.011197300241241562 + 0.021853638099321848j,
        ],
        [
            -0.01565227194902255 + 0.02289631380317063j,
            -0.010585028916259502 + 0.0008175736579655417j,
            -0.016307774327155695 + 0.027310306257580956j,
            -0.02547290396198654 - 0.013662254050241638j,
            -0.015610089840426803 + 0.022473961168034913j,
        ],
        [
            -0.0090189657144841 + 0.012557199979875446j,
            -0.01842909502111314 + 0.030862839954717714j,
            -0.003144887663914662 - 0.03142102841372029j,
            -0.017167846023984527 - 0.029801258440375466j,
            -0.00986882859033198 + 0.011381366492891239j,
        ],
        [
            -0.021558749508282765 - 0.014480490479677317j,
            -0.027838382932518536 - 0.014930965882003704j,
            -0.016602444587532073 - 0.028819791440570158j,
            -0.040690501872366926 - 0.06596298380003551j,
            -0.028162896783047266 + 0.0024728170988614065j,
        ],
        [
            -0.012781615260575384 + 0.024945726935193063j,
            -0.016051484310433087 + 0.023109440033314468j,
            -0.008979788037688476 + 0.010356067870662537j,
            -0.026498517250722525 + 0.0023266777937241134j,
            -0.007074728117953598 + 0.003313712517115923j,
        ],
    ]
    expect[1, 0, :, :] = [
        [
            -0.011258679597676765 + 0.022395921742739886j,
            -0.014099858804436134 + 0.020625426955153875j,
            -0.00718927006735804 + 0.010009695657248845j,
            -0.01777032615770584 - 0.011935898167403364j,
            -0.011197300241241562 + 0.021853638099321848j,
        ],
        [
            -0.01565227194902255 + 0.02289631380317063j,
            -0.010585028916259502 + 0.0008175736579655417j,
            -0.016307774327155695 + 0.027310306257580956j,
            -0.02547290396198654 - 0.013662254050241638j,
            -0.015610089840426803 + 0.022473961168034913j,
        ],
        [
            -0.0090189657144841 + 0.012557199979875446j,
            -0.01842909502111314 + 0.030862839954717714j,
            -0.003144887663914662 - 0.03142102841372029j,
            -0.017167846023984527 - 0.029801258440375466j,
            -0.00986882859033198 + 0.011381366492891239j,
        ],
        [
            -0.021558749508282765 - 0.014480490479677317j,
            -0.027838382932518536 - 0.014930965882003704j,
            -0.016602444587532073 - 0.028819791440570158j,
            -0.040690501872366926 - 0.06596298380003551j,
            -0.028162896783047266 + 0.0024728170988614065j,
        ],
        [
            -0.012781615260575384 + 0.024945726935193063j,
            -0.016051484310433087 + 0.023109440033314468j,
            -0.008979788037688476 + 0.010356067870662537j,
            -0.026498517250722525 + 0.0023266777937241134j,
            -0.007074728117953598 + 0.003313712517115923j,
        ],
    ]
    expect[1, 1, :, :] = [
        [
            0.979645694563844 + 0.03743860913339489j,
            -0.017857473673999053 - 0.005077635955552676j,
            -0.013734573631049728 - 0.0046245525883680996j,
            -0.024443469729725567 - 0.028669459846426687j,
            -0.018376164691808704 - 0.0016794296265635839j,
        ],
        [
            -0.019823605196670754 - 0.0056366904049734885j,
            0.9708849843014916 + 0.04031767197824951j,
            -0.0069598869114165084 - 0.05888988234327362j,
            -0.028356227993153422 - 0.05907316358840312j,
            -0.018795029420162816 - 0.031674557647168j,
        ],
        [
            -0.017230073084042043 - 0.005801518213746545j,
            -0.00786523130953012 - 0.06655029777295496j,
            0.9733774204079847 + 0.04866102775883646j,
            -0.024748282499230397 - 0.0621489255327363j,
            -0.016937172015851695 - 0.032042486369294125j,
        ],
        [
            -0.029654528360355002 - 0.034781449585204335j,
            -0.030989459802982112 - 0.06455884854994856j,
            -0.02393322891270318 - 0.06010212876387582j,
            0.9374235611255546 + 0.004370404478317288j,
            -0.021882022410520837 - 0.09240662320571531j,
        ],
        [
            -0.02097622301763208 - 0.0019170534755230829j,
            -0.019326481970050186 - 0.0325701946824621j,
            -0.015411374639662323 - 0.029155915837738133j,
            -0.020588831922819812 - 0.08694554817853925j,
            0.9720874797759053 + 0.05307586568278426j,
        ],
    ]

    def test(self):
        tm = acoustotreams.AcousticTMatrix.sphere(4, 12, [0.2], [(900, 686,), ()]).latticeinteraction.solve(
            self.basis.lattice, self.basis.kpar[:2]
        )
        sm = AcousticSMatrices.from_array(tm, self.basis)
        assert all(
            np.all(np.abs(sm[i, j] - self.expect[i, j]) < 1e-8)
            for i in range(2)
            for j in range(2)
        )

def _slab_and_waves():
    b = acoustotreams.ScalarPlaneWaveBasisByComp.default([[0, 0], [0.5, 0]])
    air = acoustotreams.AcousticMaterial()
    water = acoustotreams.AcousticMaterial(1000, 1500, 0)
    sm = AcousticSMatrices.slab(0.3, b, 2.0, [air, water, air])
    waves = [
        acoustotreams.plane_wave_scalar(kpar, k0=2.0, basis=b, material=air, modetype="up")
        for kpar in ([0, 0], [0.5, 0])
    ]
    return b, sm, waves


def test_illuminate_two_dimensional():
    b, sm, waves = _slab_and_waves()
    illu = acoustotreams.AcousticsArray(
        np.stack([np.asarray(w) for w in waves], -1),
        ({"basis": b, "modetype": "up", "k0": 2.0}, {}),
    )
    up, down = sm.illuminate(illu)
    single = [sm.illuminate(w) for w in waves]
    assert np.allclose(np.asarray(up), np.stack([np.asarray(s[0]) for s in single], -1))
    assert np.allclose(np.asarray(down), np.stack([np.asarray(s[1]) for s in single], -1))


def test_tr_two_dimensional():
    b, sm, waves = _slab_and_waves()
    illu = acoustotreams.AcousticsArray(
        np.asarray(waves[0])[:, None], ({"basis": b, "modetype": "up", "k0": 2.0}, {})
    )
    t2, r2 = sm.tr(illu)
    t1, r1 = sm.tr(waves[0])
    assert np.isclose(np.asarray(t2).item(), t1) and np.isclose(np.asarray(r2).item(), r1)


def test_tr_reordered_incident_basis():
    b, sm, waves = _slab_and_waves()
    rev = acoustotreams.ScalarPlaneWaveBasisByComp([(0.5, 0), (0, 0)])
    illu = acoustotreams.AcousticsArray(
        [0, 1], basis=rev, k0=2.0, material=acoustotreams.AcousticMaterial(), modetype="up"
    )
    assert np.allclose(sm.tr(illu), sm.tr(waves[0]))
    up, down = sm.illuminate(illu)
    up0, down0 = sm.illuminate(waves[0])
    assert np.allclose(np.asarray(up), np.asarray(up0))
    assert np.allclose(np.asarray(down), np.asarray(down0))


def test_tr_plane_wave_without_basis():
    b, sm, waves = _slab_and_waves()
    air = acoustotreams.AcousticMaterial()
    pw = acoustotreams.plane_wave_scalar([0, 0, air.ks(2.0)], k0=2.0, material=air, modetype="up")
    assert np.allclose(sm.tr(pw), sm.tr(waves[0]))


def test_slab_wrong_number_of_materials():
    b = acoustotreams.ScalarPlaneWaveBasisByComp.default([1, 2])
    with pytest.raises(ValueError):
        AcousticSMatrices.slab([3, 4], b, 12, [(200, 1029,), (900, 686,), (100, 343,)])


def test_slab_two_layers():
    b = acoustotreams.ScalarPlaneWaveBasisByComp.default([1, 2])
    mats = [(200, 1029,), (900, 686,), (900, 686,), (100, 343,)]
    two = AcousticSMatrices.slab([1, 2], b, 12, mats)
    one = AcousticSMatrices.slab(3, b, 12, [mats[0], mats[1], mats[3]])
    assert all(np.allclose(np.asarray(two[i, j]), np.asarray(one[i, j])) for i in range(2) for j in range(2))


def test_propagation_scalar_distance():
    b = acoustotreams.ScalarPlaneWaveBasisByComp.default([1, 2])
    assert AcousticSMatrices.propagation(3, b, 12) == AcousticSMatrices.propagation([0, 0, 3], b, 12)


def test_tr_several_columns():
    b, sm, waves = _slab_and_waves()
    t, r = sm.tr(np.eye(2))
    single = [sm.tr(np.eye(2)[:, i]) for i in range(2)]
    assert np.shape(t) == np.shape(r) == (2,)
    assert np.allclose(t, [s[0] for s in single]) and np.allclose(r, [s[1] for s in single])


def test_tr_evanescent_incidence():
    air = acoustotreams.AcousticMaterial()
    b = acoustotreams.ScalarPlaneWaveBasisByComp.default([[3.0, 0]])  # k = 2 in air
    sm = AcousticSMatrices.slab(0.3, b, 2.0, [air, acoustotreams.AcousticMaterial(1000, 1500, 0), air])
    with pytest.raises(ValueError):
        sm.tr([1])


def test_tr_lossy_background():
    # lossless slab between lossy media: T + R = 1
    b = acoustotreams.ScalarPlaneWaveBasisByComp.default([[0.2, 0]])
    lossy = acoustotreams.AcousticMaterial(1000, 1500 - 40j)
    slab = acoustotreams.AcousticMaterial(1050, 2350)
    sm = AcousticSMatrices.slab(0.3, b, 2.0, [lossy, slab, lossy])
    t, r = sm.tr([1])
    assert np.isclose(t + r, 1)


def test_tr_tunnelling():
    # evanescent in water, propagating in the slower medium above: all power is transmitted
    water = acoustotreams.AcousticMaterial(1000, 1500)
    slow = acoustotreams.AcousticMaterial(1200, 1100)
    b = acoustotreams.ScalarPlaneWaveBasisByComp.default([[1.15 * water.ks(2.0), 0]])
    sm = AcousticSMatrices.slab(1.0, b, 2.0, [water, water, slow])
    t, r = sm.tr([1])
    assert np.isclose(t, 1) and np.isclose(r, 0, atol=1e-12)


def test_poynting_avg_z_interference():
    # up- and down-going evanescent waves carry intensity together
    water = acoustotreams.AcousticMaterial(1000, 1500)
    k0 = 2.0
    b = acoustotreams.ScalarPlaneWaveBasisByComp.default([[1.5 * water.ks(k0), 0]])
    paz = np.asarray(acoustotreams.poynting_avg_z(b, k0, water))
    up, down = np.array([1.0]), np.array([0.4 + 0.5j])
    r = np.zeros(3)
    p = sum(
        np.asarray(acoustotreams.pfield(r, basis=b, k0=k0, material=water, modetype=m)).ravel() @ c
        for c, m in ((up, "up"), (down, "down"))
    )
    vz = sum(
        np.asarray(acoustotreams.vfield(r, basis=b, k0=k0, material=water, modetype=m))[..., 2, :].ravel() @ c
        for c, m in ((up, "up"), (down, "down"))
    )
    expect = 0.5 * np.real(p * np.conj(vz))
    got = np.real(
        up.conj() @ paz @ up - down.conj() @ paz @ down
        + down.conj() @ paz @ up - up.conj() @ paz @ down
    )
    assert expect > 0 and np.isclose(got, expect)


def test_tr_full_wave_vector_down():
    b, sm, waves = _slab_and_waves()
    air = acoustotreams.AcousticMaterial()
    down = acoustotreams.plane_wave_scalar([0, 0, -1], k0=2.0, material=air)
    ref = acoustotreams.plane_wave_scalar([0, 0, -1], k0=2.0, material=air, modetype="down")
    assert np.allclose(sm.tr(down), sm.tr(ref))


def test_illuminate_contradicting_modetype():
    b, sm, waves = _slab_and_waves()
    air = acoustotreams.AcousticMaterial()
    wrong = acoustotreams.plane_wave_scalar([0, 0, -1], k0=2.0, material=air, modetype="up")
    with pytest.raises(ValueError):
        sm.illuminate(wrong)


def test_tr_oblique_full_wave_vector():
    air = acoustotreams.AcousticMaterial()
    water = acoustotreams.AcousticMaterial(1000, 1500, 0)
    k = air.ks(2.0)
    for theta in np.linspace(0.1, 1.3, 7):
        b = acoustotreams.ScalarPlaneWaveBasisByComp.default([[k * np.sin(theta), 0]])
        sm = AcousticSMatrices.slab(0.3, b, 2.0, [air, water, air])
        pw = acoustotreams.plane_wave_scalar([np.sin(theta), 0, np.cos(theta)], k0=2.0, material=air)
        ref = acoustotreams.plane_wave_scalar([k * np.sin(theta), 0], k0=2.0, basis=b, material=air)
        assert np.allclose(sm.tr(pw), sm.tr(ref))