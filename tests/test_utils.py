from bn256.constants import P
from bn256.utils import bits_of, mod_inverse


class TestUtils:
    def test_bits_of(self):  # ✅
        k = 32498273234
        bits = bits_of(k)
        expected_bits = [1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 0, 1,
                         0, 0, 1, 0]
        assert bits == expected_bits

    def test_mod_inverse(self):
        k = 32498273234
        r = mod_inverse(k, P)
        assert r == 5113278667736460357814589262896754087238737747850571709981590827357930058526
