import numpy
import scipy
import checkerboard



def make_B(S, l, sigma):
    B = numpy.identity(S.N)
    B = right_Bmult(S, B, l, sigma)
    return B

def make_Bexact(S, l, sigma):
    B = numpy.copy(S.expT)
    B = numpy.diag(sigma * S.alpha * S.h[:, l]) @ B
    B = numpy.exp(S.dtau * S.mu) * B
    return B

def make_Bckbexact(S, l, sigma):
    B = numpy.copy(S.expTckb)
    B = numpy.diag(sigma * S.alpha * S.h[:, l]) @ B
    B = numpy.exp(S.dtau * S.mu) * B
    return B



def right_Bmult(S, A, l, sigma):
    M = A @ numpy.diag(sigma * S.aldmu * S.h[:, l])
    M = S.ckb.right_mult(M)
    return M

def right_Bmult_exact(S, A, l, sigma):
    B = make_Bexact(S, l, sigma)
    return A @ B

def right_Bmult_ckbexact(S, A, l, sigma):
    B = make_Bckbexact(S, l, sigma)
    return A @ B



def left_Bmult(S, A, l, sigma):
    M = S.ckb.left_mult(A)
    M = numpy.diag(sigma * S.aldmu * S.h[:, l]) @ M
    return M

def left_Bmult_exact(S, A, l, sigma):
    B = make_Bexact(S, l, sigma)
    return B @ A

def left_Bmult_ckbexact(S, A, l, sigma):
    B = make_Bckbexact(S, l, sigma)
    return B @ A



def right_Binvmult(S, A, l, sigma):
    M = S.ckbinv.rightmult(A)
    M = M @ numpy.diag(sigma * S.aldmuinv * S.h[:, l])
    return M

def right_Binvmult_exact(S, A, l, sigma):
    B = make_Bexact(S, l, sigma)
    return A @ scipy.linalg.inv(B)

def right_Binvmult_ckbexact(S, A, l, sigma):
    B = make_Bckbexact(S, l, sigma)
    return A @ scipy.linalg.inv(B)



def left_Binvmult(S, A, l, sigma):
    M = numpy.diag(sigma * S.aldmuinv * S.h[:, l]) @ A
    M = S.ckbinv.left_mult(M)
    return M

def left_Binvmult_exact(S, A, l, sigma):
    B = make_Bexact(S, l, sigma)
    return scipy.linalg.inv(B) @ A

def left_Binvmult_ckbexact(S, A, l, sigma):
    B = make_Bckbexact(S, l, sigma)
    return scipy.linalg.inv(B) @ A


