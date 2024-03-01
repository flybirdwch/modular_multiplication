from utils import *
from math import *

#2w + 2w = 2w+1
def machine_two_digit_add(x, y, bits_in_digit):
    """
    :param x: 2 digits number
    :param y: 2 digits number
    :param bits_in_digit: bits in each digit
    :return: x+y as 3 digit number (with carry)
    """
    X = digits_to_num(x, bits_in_digit)
    Y = digits_to_num(y, bits_in_digit)
    RES = X+Y
    res = num_to_digits(RES, bits_in_digit, 3)
    return res

# w * w = 2w
def machine_multiply(x, y, bits_in_digit):
    return num_to_digits(x*y, bits_in_digit, 2)

#多精度全乘：计算全三角：按对角线相乘后再相加，需要进位的时候进位
def mp_full_multiply(A, B, bits_in_digit, num_digits):
    if num_digits >= 2**bits_in_digit:
        raise("MP multiplier works only if: num_digits < 2**bits_in_digit")
    c = [0]*(num_digits*2 + 1)
    for l in range(num_digits*2 - 2 + 1):
        i_min = max(0, l - (num_digits - 1))
        i_max = min(l, num_digits - 1) + 1  # + 1 for inclusive
        #print('l,i_min,i_max:',l,i_min,i_max)
        for i in range(i_min, i_max):
            mult_res = machine_multiply(A[i], B[l-i], bits_in_digit)
            add_res = machine_two_digit_add(mult_res, [c[l], c[l+1]], bits_in_digit)
            [c[l], c[l+1]] = [add_res[0], add_res[1]]
            c[l+2] = c[l+2] + add_res[2] #进位
    return c[:num_digits*2]

#多精度lsb：计算上三角+对角线
def mp_lsb_multiply(A, B, bits_in_digit, num_digits, return_extra_digit=False):
    if num_digits >= 2**bits_in_digit:
        raise("MP multiplier works only if: num_digits >= 2**bits_in_digit")
    c = [0]*(num_digits*2 + 1)
    for l in range(num_digits):  # num_digits
        i_min = max(0, l - (num_digits - 1))
        i_max = min(l, num_digits - 1) + 1  # + 1 for inclusive
        #print('l,i_min,i_max:',l,i_min,i_max)
        for i in range(i_min, i_max):
            mult_res = machine_multiply(A[i], B[l-i], bits_in_digit)
            add_res = machine_two_digit_add(mult_res, [c[l], c[l+1]], bits_in_digit)
            [c[l], c[l+1]] = [add_res[0], add_res[1]]
            c[l+2] = c[l+2] + add_res[2]
    if return_extra_digit:
        return c[:num_digits+1]
    else:
        return c[:num_digits]

#计算对角线的下一条对角线
def mp_lsb_extra_diagonal(A, B, bits_in_digit, num_digits):
    if num_digits >= 2**bits_in_digit:
        raise("MP multiplier works only if: num_digits < 2**bits_in_digit")
    c = [0]*(num_digits*2 + 1)
    l = num_digits
    i_min = max(0, l - (num_digits - 1))
    i_max = min(l, num_digits - 1) + 1  # + 1 for inclusive
    print('l,i_min,i_max:',l,i_min,i_max)
    for i in range(i_min, i_max):
        mult_res = machine_multiply(A[i], B[l-i], bits_in_digit)
        add_res = machine_two_digit_add(mult_res, [c[l], c[l+1]], bits_in_digit)
        [c[l], c[l+1]] = [add_res[0], add_res[1]]
        c[l+2] = c[l+2] + add_res[2]

    return c[l]

#多精度msb：计算下三角+对角线
def mp_msb_multiply(A, B, bits_in_digit, num_digits):
    """
    Returns [A*B]_msb + e
    e is in [0, num_digits]
    """
    if num_digits >= 2**bits_in_digit:
        raise("MP multiplier works only if: num_digits < 2**bits_in_digit")
    c = [0]*(num_digits*2 + 1)
    for l in range(num_digits-1, num_digits*2 - 2 + 1):  # num_digits
        i_min = l - (num_digits - 1)
        i_max = num_digits - 1 + 1  # + 1 for inclusive
        #print('l,i_min,i_max:',l,i_min,i_max)
        for i in range(i_min, i_max):
            mult_res = machine_multiply(A[i], B[l-i], bits_in_digit)
            add_res = machine_two_digit_add(mult_res, [c[l], c[l+1]], bits_in_digit)
            [c[l], c[l+1]] = [add_res[0], add_res[1]]
            c[l+2] = c[l+2] + add_res[2]
    return c[num_digits:2*num_digits]

#+1操作
def mp_adder(A, B, bits_in_digit, num_digits):  # returns num_digits + 1 number
    C = [0] * (num_digits + 1)
    carry = 0
    for i in range(num_digits):
        carry, res = divmod(A[i] + B[i] + carry, 2**bits_in_digit)
        C[i] = res
    C[num_digits] = carry
    return C

#移位操作：将A向direction方向移动shift位
def mp_shifter(A, shift, bits_in_digit, num_digits, direction):
    a = digits_to_num(A, bits_in_digit)
    if direction == 'left':
        a_shift = a << shift
    elif direction == 'right':
        a_shift = a >> shift
    A_SHIFT = num_to_digits(a_shift, bits_in_digit, num_digits)
    return A_SHIFT

#用补码计算减法
def mp_subtract(A, B, w, n, k, z):
    a = digits_to_num(A, w)
    b = digits_to_num(B, w)
    minus_ls_plus_1 = (~b + 1) % 2 ** (n + ceil(log2(4 + (k / 2 ** z))))
    res = (a + minus_ls_plus_1) % 2 ** (n + ceil(log2(4 + (k / 2 ** z))))
    RES = num_to_digits(res, w, k+1)
    return RES

#循环减去多余的S
def mp_subtract_red(R, S, w, k):
    r = digits_to_num(R, w)
    num_red = 0
    s = digits_to_num(S, w)
    while r > s:
        r = r - s
        num_red += 1
    R = num_to_digits(r, w, k)
    return R, num_red

