from math import *
from primitive_ops import *
from utils import *
from random import *


def domb_barett_mp_redc(A, B, S, M, w, k, n):
    """
    A: input number. represented by k digits.
    B: input number. represented by k digits
    S: mod field. represented by k digits
    M: [2**(wk - n) * floor(2<<2n / s)][k-1:0] 1/S Barett approximation, represented by k bits
    w: bits in each digit
    k: number of digits required to represent S
    n: ceil(log2(S))

    outputs:
        R = A * B mod S
    """

    # Full multiply and break into LSB, MSB parts
    AB = mp_full_multiply(A, B, w, k)

    # AB msb extraction (+ shift)：从全乘取msb：先左移，再截取
    wk = w * k
    z = (wk-n)
    AB_shift = mp_shifter(AB, 2 * z, w, k * 2, 'left')
    AB_msb = AB_shift[k:2 * k]

    # L estimation
    L = mp_msb_multiply(AB_msb, M, w, k)  # calculate l estimator (MSB multiply)
    L = mp_adder(L, AB_msb, w, k)[:k]  # Add another AB_msb because m[n] = 1
    L = mp_shifter(L, z, w, k, 'right')

    # LS calculation
    LS = mp_lsb_multiply(L, S, w, k, return_extra_digit=True)
    # If needed, calculate extra diagonal.
    if z < log2(4 + k/(2**z)):
        lsb_mult_carry_extra = LS[k]
        lsb_mult_extra = mp_lsb_extra_diagonal(L, S, w, k)
        LS[k] = (lsb_mult_carry_extra + lsb_mult_extra)
    else:
        LS = LS[:k]  # remove extra digit from lsb mult

    # adders and sub, not in multiprecision.
    if z < log2(4 + k / (2 ** z)):
        AB_lsb = AB[:k + 1]
    else:
        AB_lsb = AB[:k]

    R = mp_subtract(AB_lsb, LS, w, n, k, z)
    R, num_red = mp_subtract_red(R, S, w, k)

    assert num_red < (4 + k / 2 ** z)  # check subtract redundant didn't pass bound

    return R


def domb_barrett_mp_redc_wrapper(s, a, b, w):
    n = len(bin(s)[2:]) #bin返回值开头为0x
    k = ceil(n / w)
    z = k * w - n
    m, _ = divmod(2 ** (2 * n + z), s)  # prime approximation, n + 1 bits
    A = num_to_digits(a, w, k)
    B = num_to_digits(b, w, k)
    M = num_to_digits(m, w, k + 1)[:k]
    S = num_to_digits(s, w, k)
    RES = domb_barett_mp_redc(A, B, S, M, w=w, k=k, n=n)
    res = digits_to_num(RES, w)
    return res

if __name__ == '__main__':
    #BLS12-377
    n = 377
    s = 258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458177
    w = 32 #64位的机器
    #a = randint(0, s - 1)
    #b = randint(0, s - 1)
    a = 258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458177 - 1
    b = 258664426012969094010652733694893533536393512754914660539884262666720468348340822774968888139573360124440321458177 - 1
    res = domb_barrett_mp_redc_wrapper(s, a, b, w)
    print("res", res)