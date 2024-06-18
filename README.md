# Bn256 of Python

This package is a fully tested python edition of Golang crypto/bn256.

## How to install
```shell
pip install bn256
```

## How to Use

```python
from bn256 import G1, G2, GT

k1, g1 = G1.random_g1()
k2, g2 = G2.random_g2()

gt = GT.pair(g1, g2)

print(gt.marshal().hex())
```

## Refer

- [crypto/bn256](https://pkg.go.dev/golang.org/x/crypto/bn256)
- [randombit/pairings.py](https://github.com/randombit/pairings.py/blob/master/bn256.py)