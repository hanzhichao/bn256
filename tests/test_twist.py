from bn256.gfp2 import Gfp2
from bn256.twist import TwistPoint, TWIST_G


class TestTwistPoint:

    def test_twist_point_set_infinity(self):  # ✅
        c = TwistPoint.zero().set_infinity()
        assert c.is_infinity()

    def test_twist_point_is_on_curve(self):  # ✅
        assert TWIST_G.is_on_curve() is True

    def test_twist_point_add(self):  # ✅
        c = TWIST_G + TWIST_G

        assert c.x.x == -33574719689893648050868370973934787128990408670393507348137512827186327608632
        assert c.x.y == 117509279024775130555523083373412303077470815635121613301605812856810281381884

        assert c.y.x == -41415365205126244681030222938224966145278956034555642246947848162722551516277
        assert c.y.y == -134867312344436760833210833783476190768747672672076947546926704610096902680873

        assert c.z.x == 8164735751726867362664406806290871136633702655186802416211482152428240187062
        assert c.z.y == 16991307846246862835209946494978544876836381174527200297540561298613916203860

    def test_twist_point_double(self):  # ✅
        c = TWIST_G.double()

        assert c.x.x == -33574719689893648050868370973934787128990408670393507348137512827186327608632
        assert c.x.y == 117509279024775130555523083373412303077470815635121613301605812856810281381884

        assert c.y.x == -41415365205126244681030222938224966145278956034555642246947848162722551516277
        assert c.y.y == -134867312344436760833210833783476190768747672672076947546926704610096902680873

        assert c.z.x == 8164735751726867362664406806290871136633702655186802416211482152428240187062
        assert c.z.y == 16991307846246862835209946494978544876836381174527200297540561298613916203860

    def test_twist_point_neg(self):  # ✅
        c1 = TWIST_G.negative()

        assert c1.x.x == 11559732032986387107991004021392285783925812861821192530917403151452391805634
        assert c1.x.y == 10857046999023057135944570762232829481370756359578518086990519993285655852781

        assert c1.y.x == -4082367875863433681332203403145435568316851327593401208105741076214120093531
        assert c1.y.y == -8495653923123431417604973247489272438418190587263600148770280649306958101930

        assert c1.z.x == 0
        assert c1.z.y == 1

        c2 = -TWIST_G
        assert c2.x.x == 11559732032986387107991004021392285783925812861821192530917403151452391805634
        assert c2.x.y == 10857046999023057135944570762232829481370756359578518086990519993285655852781

        assert c2.y.x == -4082367875863433681332203403145435568316851327593401208105741076214120093531
        assert c2.y.y == -8495653923123431417604973247489272438418190587263600148770280649306958101930

        assert c2.z.x == 0
        assert c2.z.y == 1

    def test_twist_point_mul_scalar(self):
        k = 32498273234
        c = TWIST_G.mul_scalar(k)
        assert c.x.x == 79342498918014555057659957993707594198152816625318490795384086866231378483238
        assert c.x.y == 75064881631888316299786248174992043733548762922009969041763335514716661580046

        assert c.y.x == -80971837679158956612470671613901869115198702943314042635428128476153813680368
        assert c.y.y == -87766621548441252636986790424712567592503460630785843968260989606226302241177

        assert c.z.x == 25982220755985358399738943490213691755613536187583364732600184316426166927358
        assert c.z.y == 28839747431195664757690418033918501226209980182353693445864133946636662806562

    def test_twist_point_mul_zero(self):
        k = 0
        c = TWIST_G.mul_scalar(k)
        assert c.x.x == 0
        assert c.x.y == 0

        assert c.y.x == 0
        assert c.y.y == 0

        assert c.z.x == 0
        assert c.z.y == 0

    def test_twist_point_make_affine(self):
        a = TwistPoint(
            Gfp2(79342498918014555057659957993707594198152816625318490795384086866231378483238,
                 75064881631888316299786248174992043733548762922009969041763335514716661580046),
            Gfp2(-80971837679158956612470671613901869115198702943314042635428128476153813680368,
                 -87766621548441252636986790424712567592503460630785843968260989606226302241177),
            Gfp2(25982220755985358399738943490213691755613536187583364732600184316426166927358,
                 28839747431195664757690418033918501226209980182353693445864133946636662806562),
        )
        c = a.copy().make_affine()
        assert c == a