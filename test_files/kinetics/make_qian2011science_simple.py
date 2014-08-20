
# write python script to make circuit.  concentrations aren't given

# use same names and numbers as Qian & Winfree, Science, 2011 supp info figure S10.  also figure S4.

# nodes are numbered by integers, preferably.  node 0 is reserved for fuel strands.

# a simplified model is used: clamps aren't modeled, and neither are extended toeholds for thresholds.

# in this directory, for the smaller circuit of figure S4:
#
# python make_qian2011science_simple.py
# more qian_andor.pil
# ../../enumerator.py qian_andor.pil -c 
# more qian_andor-enum.pil
#
# 19 reactions in the end.  painfully verified by hand.

# for the full circuit, in all its glory:
#
# python make_qian2011science_simple.py
# more qian_sqrt.pil
# ../../enumerator.py qian_sqrt.pil -c --max-complex-count 10000 --max-reaction-count 50000
# more qian_sqrt-enum.pil
#
# takes a while, but it should complete.  156 reactions after condensing.  how to check if they're right?
# (218 resting complexes, 7192 transient complexes, 14540 detailed reactions... thank you, universal toehold!)

def wire(Nfrom, Nto):
    i=str(Nto)
    j=str(Nfrom)
    return "Wire_"+i+"_to_"+j+" = S"+i+" T S"+j

def threshold(Nfrom, Nnode):
    i=str(Nnode)
    j=str(Nfrom)
    return "Threshold_"+i+" = S"+i+"( + T* )"

def gate_output(Nnode, Nto):
    i=str(Nnode)
    j=str(Nto)
    return "Gate_"+i+"_to_"+j+" = S"+j+" T( S"+i+"( + T* ))"

def reporter(name, Nnode):
    i=str(Nnode)
    return name+"_Reporter_"+i+" = S"+i+"( + T* )"

def input(name, Nfrom, Nto):
    i=str(Nto)
    j=str(Nfrom)
    return name+"_Wire_"+i+"_to_"+j+" = S"+i+" T S"+j

def fanout_gate(Nfrom, Nnode, Ntos):
    fuels = [ threshold(Nfrom,Nnode), wire(Nnode,0) ]
    for N in Ntos:
        fuels.append(gate_output(Nnode,N))
    return fuels

def ANDOR_gate(Nsum, Namp, Ntos):
    fuels = [ gate_output(Nsum,Namp), threshold(Nsum,Namp), wire(Namp,0) ]
    for N in Ntos:
        fuels.append(gate_output(Namp,N))
    return fuels

def AND_gate(Nsum, Namp, Ntos):
    return ANDOR_gate(Nsum, Namp, Ntos)

def OR_gate(Nsum, Namp, Ntos):
    return ANDOR_gate(Nsum, Namp, Ntos)

def andorcascade():
    fuels = [ input("x1",8,4), input("x2",9,4), input("x3",3,2) ]
    fuels += AND_gate(4,1,[2])
    fuels += AND_gate(2,5,[6])
    fuels.append( reporter("y_ROX",6) )
    return fuels

fuels=andorcascade()

pil_file = open("qian_andor.pil", "w")
pil_file.write("length T = 5\n")
for complex in fuels:
    pil_file.write(complex+"\n")
pil_file.close()


def squarerootcircuit():
    fuels = [ input("x01",45,42), input("x11",46,41), input("x02",47,42), input("x12",48,41), \
              input("x03",49,33), input("x13",50,35), input("x04",51,37), input("x14",52,38) ]
    fuels += fanout_gate(49,33,[34,43,26])
    fuels += fanout_gate(50,35,[36,44,20])
    fuels += fanout_gate(51,37,[34,44,26])
    fuels += fanout_gate(52,38,[36,43,20])
    fuels += OR_gate(41,28,[34,40])
    fuels += AND_gate(42,29,[36,39])
    fuels += OR_gate(43,30,[40])
    fuels += AND_gate(44,31,[39])
    fuels += OR_gate(20,8,[25])
    fuels += AND_gate(26,13,[24])
    fuels += OR_gate(34,18,[53])
    fuels += AND_gate(36,21,[10])
    fuels += OR_gate(39,22,[53])
    fuels += AND_gate(40,27,[10])
    fuels += OR_gate(10,1,[23])
    fuels += AND_gate(53,5,[6])
    fuels.append( reporter("y01_ROX",6) )
    fuels.append( reporter("y11_FAM",23) )
    fuels.append( reporter("y02_TYE563",24) )
    fuels.append( reporter("y12_TYE665",25) )
    return fuels

fuels=squarerootcircuit()

pil_file = open("qian_sqrt.pil", "w")
pil_file.write("length T = 5\n")
for complex in fuels:
    pil_file.write(complex+"\n")
pil_file.close()


