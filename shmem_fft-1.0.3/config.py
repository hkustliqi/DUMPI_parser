#!/usr/bin/python

#                   UNCLASSIFIED//FOR OFFICIAL USE ONLY                       #

import sys, math

if len(sys.argv) != 5:
    print "config.py node_count cores_per_node GiBytes_per_node fraction_of_memory"
    sys.exit(1)

nodes          = int(sys.argv[1])
cores_per_node = int(sys.argv[2])
bytes_per_node = float(sys.argv[3]) * (1 << 30)
percent        = float(sys.argv[4])

max_memory     = bytes_per_node * nodes * percent / 2   # benchmark requires two large arrays

def LK( f, W ):
    ( w, z ) = ( i / 2 for i in f )

    L  = 2**w * 3**(f[1] - z)
    K  = 2**(f[0] - w) * 3**z

    lw = (L + W - 1) / W
    kw = (K + W - 1) / W

    l  = W * lw
    k  = W * kw

    return ( L, K, l, k )

def footprint( nodes, cores_per_node, f ):
    ( L, K, l, k ) = LK( f, nodes * cores_per_node )
    n = l * k

    return n * 2 * 4                  # complex floats are 8 bytes

def sheap( cores_per_node, memory  ):
    return "export XT_SYMMETRIC_HEAP_SIZE=%.0fm" % math.ceil(memory / cores_per_node / (1 << 20) )

def runline( nodes, cores_per_node, f ):
    ( w, z ) = ( i / 2 for i in f )
    ( L, K ) = LK( f, nodes * cores_per_node )[0:2]

    return "aprun -N %u -n %u ./sfft %u %u 4 1" % ( cores_per_node, nodes * cores_per_node, L, K )

f = ( 0, 0 )

while footprint( nodes, cores_per_node, f ) <= max_memory:
    l = f
    f = ( f[0] + 1,  f[1] )
f = l

print sheap( cores_per_node, percent * bytes_per_node )

print runline( nodes, cores_per_node, f )

while f[0] > 0:
    f = ( f[0], f[1] + 1 )

    while f[0] > 0 and footprint( nodes, cores_per_node, f ) > max_memory:
        f = ( f[0] - 1, f[1] )

    if f[0] >= 0:
        print runline( nodes, cores_per_node, f )


#                     UNCLASSIFIED//FOR OFFICIAL USE ONLY                     #
