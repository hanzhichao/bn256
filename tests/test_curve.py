from bn256.curve import CurvePoint
from bn256.utils import bits_of


class TestCurvePoint:
    a = CurvePoint(x=1, y=-2, z=1)
    b = CurvePoint(x=79885311972705142798326482969936249219924770158001168883491309517089224520499,
                   y=-117978995929271306268700300631116017955026219006325540957535965864796039810533,
                   z=33296282955968814767393647671175158596516707523800561594047204110963892554884)
    k = 32498273234

    def test_curve_point_init(self):
        assert self.b.x == 79885311972705142798326482969936249219924770158001168883491309517089224520499
        assert self.b.y == -117978995929271306268700300631116017955026219006325540957535965864796039810533
        assert self.b.z == 33296282955968814767393647671175158596516707523800561594047204110963892554884

    def test_curve_point_add(self):
        c1 = self.a + self.b
        assert c1.x == -8030019297004030839387309015943663447814033459498803802815950552962657336198
        assert c1.y == 17369015046471995974106459814434955140906951195137422436589153304383829678254
        assert c1.z == 17258309029904047582215572897898954019212799630461057332267253245789321192076
        assert c1.t == 0

        c2 = self.b.add(self.a)
        assert c2.x == -29918262168843306061633714761200938536510344616796627465504988447607883544781
        assert c2.y == -17369015046471995974106459814434955140906951195137422436589153304383829678254
        assert c2.z == 4629933841935227640030832847358321069483511526836766330421784648855905016507
        assert c2.t == 0

        assert c2 == c1

    def test_curve_point_make_affine(self):  # ✅
        c = self.b.copy()
        c = c.make_affine()
        assert c.x == 7483470414448436599363905724866355193253920941172561288805573139617879816370
        assert c.y == 128747707450769087512959846171490976965044225949551331963869426418068834555
        assert c.z == 1
        assert c.t == 1

    def test_curve_point_string(self):  # ✅
        c = self.b.copy()
        assert c.string() == ('(7483470414448436599363905724866355193253920941172561288805573139617879816370,'
                              '128747707450769087512959846171490976965044225949551331963869426418068834555)')

    def test_curve_point_repr(self):  # ✅
        c = self.b.copy()
        assert repr(c) == ('<CurvePoint (79885311972705142798326482969936249219924770158001168883491309517089224520499,'
                           '-117978995929271306268700300631116017955026219006325540957535965864796039810533,'
                           '33296282955968814767393647671175158596516707523800561594047204110963892554884)>')

    def test_curve_point_copy(self):  # ✅
        c = self.b.copy()
        assert c == self.b

    def test_curve_point_double(self):  # ✅
        c = self.b.double()

        assert c.x == 935411005489017982444253869730368550420232897228501539700537887390341768186
        assert c.y == -63329278205982424097876106502588387707621067484724644510142425134433553997206
        assert c.z == 20126853473059445258968919811389535708805351903341431478509918543146689543950
        assert c.t == 0

    def test_curve_point_zero_double(self):  # ✅
        c = CurvePoint.zero().double()

        assert c.x == 0
        assert c.y == 0
        assert c.z == 0
        assert c.t == 0

    def test_curve_point_mul_scalar(self):  # ✅
        c = self.a.mul_scalar(self.k)

        assert c.x == 79885311972705142798326482969936249219924770158001168883491309517089224520499
        assert c.y == -117978995929271306268700300631116017955026219006325540957535965864796039810533
        assert c.z == 33296282955968814767393647671175158596516707523800561594047204110963892554884
        assert c.t == 0

    def test_curve_point_mul_zero(self):  # ✅

        c = self.a.mul_scalar(0)

        assert c.x == 0
        assert c.y == 0
        assert c.z == 0
        assert c.t == 0

    def test_curve_point_is_on_curve(self):  # ✅
        assert self.a.is_on_curve() is True
        assert self.b.is_on_curve() is False

    def test_curve_point_is_infinite(self):  # ✅
        c = self.b.copy()
        assert c.is_infinity() is False
        c.set_infinity()
        assert c.is_infinity() is True

    def test_curve_point_set(self):  # ✅
        c = CurvePoint.zero().set(79885311972705142798326482969936249219924770158001168883491309517089224520499,
                                  -117978995929271306268700300631116017955026219006325540957535965864796039810533,
                                  33296282955968814767393647671175158596516707523800561594047204110963892554884)
        assert c == self.b

    def test_curve_point_negative(self):  # ✅
        c = -self.b
        assert c.x == 79885311972705142798326482969936249219924770158001168883491309517089224520499
        assert c.y == 117978995929271306268700300631116017955026219006325540957535965864796039810533
        assert c.z == 33296282955968814767393647671175158596516707523800561594047204110963892554884
        assert c.t == 0

        c2 = self.b.negative()
        assert c == c2

    def test_curve_point_double_2(self):  # ✅
        sum = CurvePoint(-43776485743495046546487551583707039415829214237436370882991831551060397914766,
                         -9467921253766007622712243895957299289870386820263329380995628705922624299036,
                         350834808772941454898449895915520)
        c = sum.double()
        assert c.x == 10790957796095752226823453676166352787855536353631975827840056415610578012430
        assert c.y == -14554534243840831968271509262327092318269985762466049132223060768894044758975
        assert c.z == 37546741787101857219979687028463595807101726323092929554000333841513965935888
        assert c.t == 0
