# from mpi4py import MPI
# import os

# print(os.cpu_count())

# comm = MPI.COMM_WORLD
# rank = comm.Get_rank()

# print("This is rank {}".format(rank))

def F(a):
    return(a,a+1)

(a,b) = F(0)

print(a,b)