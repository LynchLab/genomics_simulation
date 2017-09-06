#!/software/bin/python-2.7.9

import sys
import numpy

mn=sys.argv[1]
mx=sys.argv[2]

print numpy.random.uniform(mn, mx)
