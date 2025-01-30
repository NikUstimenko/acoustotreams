import numpy as np
import pytest

import acoustotreams


def isclose(a, b, rel_tol=1e-09, abs_tol=0.0):
    return abs(a - b) <= max(rel_tol * max(abs(a), abs(b)), abs_tol)


class TestSSWB:
    def test_init_empty(self):
        b = acoustotreams.ScalarSphericalWaveBasis([])
        assert b.l.size == 0 and b.m.size == 0 and b.pidx.size == 0

    def test_init_numpy(self):
        b = acoustotreams.ScalarSphericalWaveBasis(np.array([[0, 0]]), [0, 1, 0])
        assert (
            np.all(b.l == [0])
            and np.all(b.m == [0])
            and np.all(b.pidx == [0])
        )

    def test_init_duplicate(self):
        b = acoustotreams.ScalarSphericalWaveBasis([[0, 0, 0], [0, 0, 0]])
        assert (
            np.all(b.l == [0])
            and np.all(b.m == [0])
            and np.all(b.pidx == [0])
        )

    def test_init_invalid_shape(self):
        with pytest.raises(ValueError):
            acoustotreams.ScalarSphericalWaveBasis([[0], [0, 0]])

    def test_init_invalid_positions(self):
        with pytest.raises(ValueError):
            acoustotreams.ScalarSphericalWaveBasis([[0, 0]], [1, 2])

    def test_init_non_int_value(self):
        with pytest.raises(ValueError):
            acoustotreams.ScalarSphericalWaveBasis([[0.1, 0]])

    def test_init_negative_l(self):
        with pytest.raises(ValueError):
            acoustotreams.ScalarSphericalWaveBasis([[-1, 0]])

    def test_init_non_too_large_m(self):
        with pytest.raises(ValueError):
            acoustotreams.ScalarSphericalWaveBasis([[1, -2]])

    def test_init_unspecified_positions(self):
        with pytest.raises(ValueError):
            acoustotreams.ScalarSphericalWaveBasis([[1, 1, 0]])

    def test_property_positions(self):
        a = np.array([[1, 2, 3]])
        b = acoustotreams.ScalarSphericalWaveBasis([[0, 0]], a)
        assert (a == b.positions).all()

    def test_repr(self):
        b = acoustotreams.ScalarSphericalWaveBasis(
            [[0, 0, 0], [1, 0, 0]], [[0.0, 0, 0], [1, 0, 0]]
        )
        assert (
            repr(b)
            == """ScalarSphericalWaveBasis(
    pidx=[0 1],
    l=[0 0],
    m=[0 0],
    positions=[[0. 0. 0.], [1. 0. 0.]],
)"""
        )

    def test_getitem_plm(self):
        b = acoustotreams.ScalarSphericalWaveBasis([[0, 2, -1]])
        assert b.plm == ([0], [2], [-1])

    def test_getitem_lms(self):
        b = acoustotreams.ScalarSphericalWaveBasis([[2, -1]])
        assert b.lm == ([2], [-1])

    def test_getitem_invalid_index(self):
        b = acoustotreams.ScalarSphericalWaveBasis([[0, 0]])
        with pytest.raises(AttributeError):
            b.fail

    def test_getitem_int(self):
        b = acoustotreams.ScalarSphericalWaveBasis([[1, 0], [1, -1]])
        assert b[1] == (0, 1, -1)

    def test_getitem_tuple(self):
        b = acoustotreams.ScalarSphericalWaveBasis([[1, 0], [1, -1]])
        assert (np.array(b[()]) == ([0, 0], [1, 1], [0, -1])).all()

    def test_getitem_slice(self):
        a = acoustotreams.ScalarSphericalWaveBasis([[1, 0]])
        b = acoustotreams.ScalarSphericalWaveBasis([[1, 0], [1, 1]])
        assert a == b[:1]

    def test_default(self):
        a = acoustotreams.ScalarSphericalWaveBasis.default(2, 2, [[0, 0, 0], [1, 0, 0]])
        b = acoustotreams.ScalarSphericalWaveBasis(
            zip(
                9 * [0] + 9 * [1],
                2 * ([0] + 3 * [1] + 5 * [2]),
                2 * [0, -1, 0, 1, -2, -1, 0, 1, 2],
            ),
            [[0, 0, 0], [1, 0, 0]],
        )
        assert a == b

    def test_defaultlmax(self):
        assert acoustotreams.ScalarSphericalWaveBasis.defaultlmax(32, 2) == 3

    def test_defaultlmax_fail(self):
        with pytest.raises(ValueError):
            acoustotreams.ScalarSphericalWaveBasis.defaultlmax(2)

    def test_defaultdim(self):
        assert acoustotreams.ScalarSphericalWaveBasis.defaultlmax(3, 2) == 32

    def test_defaultdim_fail(self):
        with pytest.raises(ValueError):
            acoustotreams.ScalarSphericalWaveBasis.defaultdim(1, -1)

    def test_property_isglobal_true(self):
        b = acoustotreams.ScalarSphericalWaveBasis([[0, 0]])
        assert b.isglobal

    def test_property_isglobal_false(self):
        b = acoustotreams.ScalarSphericalWaveBasis(
            [[0, 0, 0], [1, 0, 0]], [[0, 0, 0], [1, 0, 0]]
        )
        assert not b.isglobal

    def test_from_iterable(self):
        a = acoustotreams.ScalarSphericalWaveBasis.default(1)
        b = acoustotreams.ScalarSphericalWaveBasis.default(2)
        assert a & b == a

    def test_neq(self):
        b = acoustotreams.ScalarSphericalWaveBasis.default(1)
        assert not b == []


class TestSCWB:
    def test_init_empty(self):
        b = acoustotreams.ScalarCylindricalWaveBasis([])
        assert b.kz.size == 0 and b.m.size == 0 and b.pidx.size == 0

    def test_init_numpy(self):
        b = acoustotreams.ScalarCylindricalWaveBasis(np.array([[0.5, 0]]), [0, 1, 0])
        assert (
            np.all(b.kz == [0.5])
            and np.all(b.m == [0])
            and np.all(b.pidx == [0])
        )

    def test_init_duplicate(self):
        b = acoustotreams.ScalarCylindricalWaveBasis([[0, 1, 0], [0, 1, 0]])
        assert (
            np.all(b.kz == [1])
            and np.all(b.m == [0])
            and np.all(b.pidx == [0])
        )

    def test_init_invalid_shape(self):
        with pytest.raises(ValueError):
            acoustotreams.ScalarCylindricalWaveBasis([[0], [0, 0]])

    def test_init_invalid_positions(self):
        with pytest.raises(ValueError):
            acoustotreams.ScalarCylindricalWaveBasis([[1, 0]], [1, 2])

    def test_init_non_int_value(self):
        with pytest.raises(ValueError):
            acoustotreams.ScalarCylindricalWaveBasis([[1, 0.2]])

    def test_init_unspecified_positions(self):
        with pytest.raises(ValueError):
            acoustotreams.ScalarCylindricalWaveBasis([[1, 1, 0]])

    def test_property_positions(self):
        a = np.array([[1, 2, 3]])
        b = acoustotreams.ScalarCylindricalWaveBasis([[1, 0]], a)
        assert (a == b.positions).all()

    def test_repr(self):
        b = acoustotreams.ScalarCylindricalWaveBasis(
            [[0, 1, 0], [1, 1, 0]], [[0.0, 0, 0], [1, 0, 0]]
        )
        assert (
            repr(b)
            == """ScalarCylindricalWaveBasis(
    pidx=[0 1],
    kz=[1. 1.],
    m=[0 0],
    positions=[[0. 0. 0.], [1. 0. 0.]],
)"""
        )

    def test_getitem_pzm(self):
        b = acoustotreams.ScalarCylindricalWaveBasis([[0, 2, -1]])
        assert b.pzms == ([0], [2], [-1])

    def test_getitem_zm(self):
        b = acoustotreams.ScalarCylindricalWaveBasis([[2, -1]])
        assert b.zms == ([2], [-1])

    def test_getitem_invalid_index(self):
        b = acoustotreams.ScalarCylindricalWaveBasis([[1, 0]])
        with pytest.raises(AttributeError):
            b.fail

    def test_getitem_int(self):
        b = acoustotreams.ScalarCylindricalWaveBasis([[1, 0], [1, -1]])
        assert b[0] == (0, 1.0, -1)

    def test_getitem_tuple(self):
        b = acoustotreams.ScalarCylindricalWaveBasis([[1, 0], [1, -1]])
        assert (np.array(b[()]) == ([0., 0.], [1., 1.], [0, -1.])).all()

    def test_getitem_slice(self):
        a = acoustotreams.ScalarCylindricalWaveBasis([[1, 0]])
        b = acoustotreams.ScalarCylindricalWaveBasis([[1, 0], [1, 1]])
        assert a == b[:1]

    def test_default(self):
        a = acoustotreams.ScalarCylindricalWaveBasis.default(
            [0.3, -0.2], 2, 2, [[0, 0, 0], [1, 0, 0]]
        )
        b = acoustotreams.ScalarCylindricalWaveBasis(
            zip(
                10 * [0] + 10 * [1],
                2 * (5 * [0.3] + 5 * [-0.2]),
                4 * [-2, -1, 0, 1, 2],
            ),
            [[0, 0, 0], [1, 0, 0]],
        )
        assert a == b

    def test_defaultmmax(self):
        assert acoustotreams.ScalarCylindricalWaveBasis.defaultmmax(56, 4, 2) == 3

    def test_defaultmmax_fail(self):
        with pytest.raises(ValueError):
            acoustotreams.ScalarCylindricalWaveBasis.defaultmmax(2)

    def test_defaultdim(self):
        assert acoustotreams.ScalarCylindricalWaveBasis.defaultdim(4, 3, 2) == 56

    def test_defaultdim_fail(self):
        with pytest.raises(ValueError):
            acoustotreams.ScalarCylindricalWaveBasis.defaultdim(1, -1)

    def test_property_isglobal_true(self):
        b = acoustotreams.ScalarCylindricalWaveBasis.defaultdim([[1, 0], [1, 1]])
        assert b.isglobal

    def test_property_isglobal_false(self):
        b = acoustotreams.ScalarCylindricalWaveBasis(
            [[0, 1, 0], [1, 1, 0]], [[0, 0, 0], [1, 0, 0]]
        )
        assert not b.isglobal

    def test_from_iterable(self):
        a = acoustotreams.ScalarCylindricalWaveBasis(0, 1)
        b = acoustotreams.ScalarCylindricalWaveBasis(0, 2)
        assert a & b == a

    def test_neq(self):
        b = acoustotreams.ScalarCylindricalWaveBasis(0, 1)
        assert not b == []

    def test_diffr_orders(self):
        a = acoustotreams.ScalarCylindricalWaveBasis.diffr_orders(0.1, 1, 2 * np.pi, 1.5)
        b = acoustotreams.ScalarCylindricalWaveBasis(
            zip(
                *[
                    3 * [-0.9] + 3 * [0.1] + 3 * [1.1],
                    3 * [-1, 0, 1],
                ]
            )
        )
        assert (
            a == b
            and a.lattice == acoustotreams.Lattice(2 * np.pi)
            and a.kpar == [np.nan, np.nan, 0.1]
        )


class TestSPWBUV:
    def test_init_empty(self):
        b = acoustotreams.ScalarPlaneWaveBasisByUnitVector([])
        assert b.qx.size == 0 and b.qy.size == 0 and b.qz.size == 0

    def test_init_numpy(self):
        b = acoustotreams.ScalarPlaneWaveBasisByUnitVector(np.array([[0.1, 0.2, 0.2]]))
        assert (
            np.all(b.qx == [1 / 3])
            and np.all(b.qy == [2 / 3])
            and np.all(b.qz == [2 / 3])
        )

    def test_init_duplicate(self):
        b = acoustotreams.ScalarPlaneWaveBasisByUnitVector([[0.1, 0.2, 0.2], [0.1, 0.2, 0.2]])
        assert (
            np.all(b.qx == [1 / 3])
            and np.all(b.qy == [2 / 3])
            and np.all(b.qz == [2 / 3])
        )

    def test_init_invalid_shape(self):
        with pytest.raises(ValueError):
            acoustotreams.ScalarPlaneWaveBasisByUnitVector([[0, 0], [0, 0]])

    def test_repr(self):
        b = acoustotreams.ScalarPlaneWaveBasisByUnitVector([[0.0, 1.0, 0.0], [1, 0, 0]])
        assert (
            repr(b)
            == """ScalarPlaneWaveBasisByUnitVector(
    qx=[0. 1.],
    qy=[1. 0.],
    qz=[0. 0.],
)"""
        )

    def test_getitem_xys(self):
        b = acoustotreams.ScalarPlaneWaveBasisByUnitVector([[0, 4, -3]])
        assert b.xys == ([0], [0.8], [1])

    def test_getitem_invalid_index(self):
        b = acoustotreams.ScalarPlaneWaveBasisByUnitVector([[1, 0, 0]])
        with pytest.raises(AttributeError):
            b.fail

    def test_getitem_int(self):
        b = acoustotreams.ScalarPlaneWaveBasisByUnitVector([[1, 0, 0], [0, 0, 1]])
        assert b[1] == (0., 0., 1.)

    def test_getitem_tuple(self):
        b = acoustotreams.ScalarPlaneWaveBasisByUnitVector([[1, 0, 0], [0, 0, 1]])
        assert (np.array(b[()]) == ([1, 0], [0, 0], [0, 1])).all()

    def test_getitem_slice(self):
        a = acoustotreams.ScalarPlaneWaveBasisByUnitVector([[1, 0, 0]])
        b = acoustotreams.ScalarPlaneWaveBasisByUnitVector([[1, 0, 0], [0, 0, 1]])
        assert a == b[:1]

    def test_default(self):
        a = acoustotreams.ScalarPlaneWaveBasisByUnitVector.default([0.3, -0.2, 0.1])
        b = acoustotreams.ScalarPlaneWaveBasisByUnitVector(
            [[0.3, -0.2, 0.1]]
        )
        assert a == b

    def test_property_isglobal_true(self):
        b = acoustotreams.ScalarPlaneWaveBasisByUnitVector([])
        assert b.isglobal

    def test_from_iterable(self):
        a = acoustotreams.ScalarPlaneWaveBasisByUnitVector.default([0, 0, 1])
        b = acoustotreams.ScalarPlaneWaveBasisByUnitVector.default([[0, 0, 1], [0, 1, 0]])
        assert a & b == a

    def test_bycomp(self):
        a = acoustotreams.ScalarPlaneWaveBasisByComp.default([0, 1], "yz")
        b = acoustotreams.ScalarPlaneWaveBasisByUnitVector.default([0, 0, 1])
        assert b.bycomp(1, "yz") == a


class TestSPWBC:
    def test_init_empty(self):
        b = acoustotreams.ScalarPlaneWaveBasisByComp([])
        assert b.kx.size == 0 and b.ky.size == 0 and b.kz is None

    def test_init_numpy(self):
        b = acoustotreams.ScalarPlaneWaveBasisByComp(np.array([[0.4, 0.2]]), "zx")
        assert (
            np.all(b.kz == [0.4])
            and np.all(b.kx == [0.2])
            and np.all(b.ky is None)
        )

    def test_init_duplicate(self):
        b = acoustotreams.ScalarPlaneWaveBasisByComp([[0.4, 0.2], [0.4, 0.2]])
        assert (
            np.all(b.kx == [0.4])
            and np.all(b.ky == [0.2])
            and np.all(b.kz is None)
        )

    def test_init_invalid_shape(self):
        with pytest.raises(ValueError):
            acoustotreams.ScalarPlaneWaveBasisByComp([[0], [0]])


    def test_repr(self):
        b = acoustotreams.ScalarPlaneWaveBasisByComp([[0.0, 1.0], [1, 1]], "yz")
        assert (
            repr(b)
            == """ScalarPlaneWaveBasisByComp(
    ky=[0. 1.],
    kz=[1. 1.],
)"""
        )

    def test_from_iterable(self):
        a = acoustotreams.ScalarPlaneWaveBasisByComp.default([0, 0])
        b = acoustotreams.ScalarPlaneWaveBasisByComp.default([[0, 0], [0, 1]])
        assert a & b == a

    def test_getitem_int(self):
        b = acoustotreams.ScalarPlaneWaveBasisByComp([[1, 0], [0, 1]])
        assert b[1] == (0, 1)

    def test_getitem_tuple(self):
        b = acoustotreams.ScalarPlaneWaveBasisByComp([[1, 0], [0, 1]])
        assert (np.array(b[()]) == ([1, 0], [0, 1])).all()

    def test_getitem_slice(self):
        a = acoustotreams.ScalarPlaneWaveBasisByComp([[1, 0]])
        b = acoustotreams.ScalarPlaneWaveBasisByComp([[1, 0], [0, 1]])
        assert a == b[:1]

    def test_byunitvector(self):
        a = acoustotreams.ScalarPlaneWaveBasisByUnitVector.default([0, 0, 1])
        b = acoustotreams.ScalarPlaneWaveBasisByComp.default([0, 1], "yz")
        assert b.byunitvector(1) == a

    def test_diffr_orders(self):
        lattice = acoustotreams.Lattice(2 * np.pi * np.eye(2))
        b = acoustotreams.ScalarPlaneWaveBasisByComp.diffr_orders([0, 0], lattice, 1)
        a = acoustotreams.ScalarPlaneWaveBasisByComp.default(
            [[0, 0], [0, 1], [1, 0], [-1, 0], [0, -1]]
        )
        assert a <= b and b <= a and b.lattice == lattice and b.kpar == [0, 0, np.nan]


class TestAcousticsArray:
    def test_init(self):
        b = acoustotreams.ScalarSphericalWaveBasis([[0, 0], [1, 0]])
        p = acoustotreams.AcousticsArray(np.eye(2), basis=b)
        assert (p == np.eye(2)).all() and p.basis == b

    def test_type_error(self):
        with pytest.raises(TypeError):
            acoustotreams.AcousticsArray(np.eye(2), basis="fail")

    def test_lattice(self):
        b = acoustotreams.ScalarPlaneWaveBasisByComp.diffr_orders([0, 0], np.eye(2), 4)
        p = acoustotreams.AcousticsArray([1, 2], lattice=acoustotreams.Lattice(1, "x"), basis=b)
        assert p.lattice == acoustotreams.Lattice(1, "x")

    def test_matmul(self):
        b = acoustotreams.ScalarSphericalWaveBasis([[0, 0]])
        p = acoustotreams.AcousticsArray([[1, 1], [2, 2]], basis=b)
        x = p @ [1, 2]
        assert (x == [3, 6]).all() and x.basis == b

    def test_rmatmul(self):
        b = acoustotreams.ScalarSphericalWaveBasis([[0, 0]])
        p = acoustotreams.AcousticsArray([[1, 1], [2, 2]], basis=b)
        x = [1, 2] @ p
        assert (x == [5, 5]).all() and x.basis == b

        #do others