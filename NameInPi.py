"""
Python3 program to calculate Pi using python long integers, binary
splitting and the Chudnovsky algorithm

See: http://www.craig-wood.com/nick/articles/pi-chudnovsky/ for more
info

Nick Craig-Wood <nick@craig-wood.com>
"""

import math
import gmpy2
from gmpy2 import mpz
from time import time
import argparse

def check_string(string, number):
    number = str(number)
    enum_string = enumerate(string)
    len_string = len(string)
    len_number = len(number)
    j = 0
    for i in range(len(number)):
        for j in range(len(string)):
            if number[i+j*3:i+j*3+3] == '':
                return -1
            # print(number[i+j*3:i+j*3+3])
            # print(string[j], ' =? ', str(chr(int(number[i+j*3:i+j*3+3]))))
            decoded_char = str(chr(int(number[i+j*3:i+j*3+3])))
            if decoded_char == string[j]:
                if j == len_string-1:
                    return i
            else:
                break
    return -1


def pi_chudnovsky_bs(digits):
    """
    Compute int(pi * 10**digits)

    This is done using Chudnovsky's series with binary splitting
    """
    C = 640320
    C3_OVER_24 = C**3 // 24
    def bs(a, b):
        """
        Computes the terms for binary splitting the Chudnovsky infinite series

        a(a) = +/- (13591409 + 545140134*a)
        p(a) = (6*a-5)*(2*a-1)*(6*a-1)
        b(a) = 1
        q(a) = a*a*a*C3_OVER_24

        returns P(a,b), Q(a,b) and T(a,b)
        """
        if b - a == 1:
            # Directly compute P(a,a+1), Q(a,a+1) and T(a,a+1)
            if a == 0:
                Pab = Qab = mpz(1)
            else:
                Pab = mpz((6*a-5)*(2*a-1)*(6*a-1))
                Qab = mpz(a*a*a*C3_OVER_24)
            Tab = Pab * (13591409 + 545140134*a) # a(a) * p(a)
            if a & 1:
                Tab = -Tab
        else:
            # Recursively compute P(a,b), Q(a,b) and T(a,b)
            # m is the midpoint of a and b
            m = (a + b) // 2
            # Recursively calculate P(a,m), Q(a,m) and T(a,m)
            Pam, Qam, Tam = bs(a, m)
            # Recursively calculate P(m,b), Q(m,b) and T(m,b)
            Pmb, Qmb, Tmb = bs(m, b)
            # Now combine
            Pab = Pam * Pmb
            Qab = Qam * Qmb
            Tab = Qmb * Tam + Pam * Tmb
        return Pab, Qab, Tab
    # how many terms to compute
    DIGITS_PER_TERM = math.log10(C3_OVER_24/6/2/6)
    N = int(digits/DIGITS_PER_TERM + 1)
    # Calclate P(0,N) and Q(0,N)
    P, Q, T = bs(0, N)
    one_squared = mpz(10)**(2*digits)
    sqrtC = sqrtC = gmpy2.isqrt(10005*one_squared)
    return (Q*426880*sqrtC) // T

# The last 5 digits or pi for various numbers of digits
check_digits = {
        100 : 70679,
       1000 :  1989,
      10000 : 75678,
     100000 : 24646,
    1000000 : 58151,
   10000000 : 55897,
}

if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Find the position where your word is encoded in the number pi")
    parser.add_argument("word", help="word to be found",
                        type=str)
    args = parser.parse_args()

    digits = 100
    pi = pi_chudnovsky_bs(digits)
    print(pi)
    place = check_string(args.word, pi)
    if place > 0:
        print(args.word, ' is found at ', place, ' in the pi digit!')
        exit(0)
    #raise SystemExit
    for log10_digits in range(100):
        digits = 10**log10_digits
        start =time()
        pi = pi_chudnovsky_bs(digits)
        print("chudnovsky_gmpy_mpz_bs: digits",digits,"time",time()-start)
        place = check_string(args.word, pi)
        if place > 0:
            print(args.word, ' is found at ', place, ' in the pi digit!')
            exit(0)
        # if digits in check_digits:
        #     last_five_digits = pi % 100000
        #     if check_digits[digits] == last_five_digits:
        #         print("Last 5 digits %05d OK" % last_five_digits)
        #     else:
        #         print("Last 5 digits %05d wrong should be %05d" % (last_five_digits, check_digits[digits]))
