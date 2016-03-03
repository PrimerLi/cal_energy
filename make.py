#!/usr/bin/env python 

def main():
    import os
    import sys

    os.system("gfortran -c energy.f90")
    os.system("gfortran -c fit.f90")
    os.system("gfortran -o energy.out energy.o fit.o")
    return 0

if __name__ == "__main__":
    import sys
    sys.exit(main())
