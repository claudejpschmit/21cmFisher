import sys
import numpy as np
import operator

filename_in = str(sys.argv[1])
if len(sys.argv) == 3:
    filename_out = str(sys.argv[2])
else:
    filename_out = filename_in

def checkline(line):
    if (line == "\n") or (line == ""):
        return False
    else:
        return True
    
with open(filename_in, "r") as f:
    sorted_file = sorted(filter(checkline,f))

with open(filename_out, "w+") as f:
    for line in sorted_file:
        f.write(line)
