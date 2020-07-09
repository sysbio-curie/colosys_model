# Translates a .bnet file describing a logical model into two MaBoSS files (.bnd and .cfg)

import sys

if len(sys.argv)<2:
    print "Usage: python bnet_to_bnd.py file.bnet"
    sys.exit(1)

bnetfile=sys.argv[1]

text=''
noderulelist=list()
fo = open(bnetfile,"r")
for line in fo:
    line=line.replace("\n","")
    line=line.replace("\r","")
    line=line.replace(" ","")
    node, rule=line.split(",")
    if node!="targets": # if the .bnet file contains a header "targets, factors"
        noderulelist.append((node,rule))
        text+="Node "+node+" {\n  logic = "+rule+";\n  rate_up = @logic ? $u_"+node+" : 0;\n  rate_down = @logic ? 0 : $d_"+node+";\n}\n\n"

fo.close()


## -- Generate .bnd --

bndfile=bnetfile.replace("bnet","bnd")
fo2 = open(bndfile,"w")
fo2.write(text)
fo2.close()

print("Translated "+bnetfile+" into "+bndfile)


## -- Generate .cfg --

cfgfile=bnetfile.replace("bnet","cfg")
fo3 = open(cfgfile,"w")

# for each node, add parameter lines: $u_node=1;\n$d_node=1; 
for node,rule in noderulelist:
     fo3.write("$u_"+node+"=1;\n$d_"+node+"=1;\n")

fo3.write("\n")
# add initial states for constant or input nodes:
# target, 1 -> target.istate=1;
# target, 0 -> target.is_state=0;
# target, target -> target.is_state=0;
for node,rule in noderulelist:
     if rule=="1":
         fo3.write(node+".istate=1;\n")
     elif rule=="0":
         fo3.write(node+".istate=0;\n")
     elif rule=="("+node+")" or rule==node:
         fo3.write(node+".istate=0;\n")

fo3.write("\n")
# for each node add parameter line: node.is_internal=TRUE;
for node,rule in noderulelist:
     fo3.write(node+".is_internal=TRUE;\n")

# simulations parameters
fo3.write("\ndiscrete_time=0;\nuse_physrandgen=FALSE;\nseed_pseudorandom=100;\nsample_count=200;\nmax_time=30;\ntime_tick=0.01;\nstatdist_traj_count=100;\nstatdist_cluster_threshold=0.9;\nthread_count = 4;\n")

fo3.close()

print("Translated "+bnetfile+" into "+cfgfile)

