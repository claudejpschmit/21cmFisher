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

content_list = []
for s in sorted_file:
    content_list.append(s.split())

first_elem_list = []
for l in content_list:
    first_elem_list.append(l[0])

dictionary = {l[0]: [l[1],l[2]] for l in content_list}
sorted_fe = sorted(first_elem_list, key = int)

out_file = []
for elem in sorted_fe:
    res = elem + " " + dictionary[elem][0] + " " + dictionary[elem][1] + "\n"
    out_file.append(res)

with open(filename_out, "w+") as f:
    for line in out_file:
        f.write(line)
