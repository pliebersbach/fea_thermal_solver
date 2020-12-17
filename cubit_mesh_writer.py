#!python
import os

num_nodes = cubit.get_node_count()
num_elements = cubit.get_quad_count()
num_edges = cubit.get_edge_count()

num_nodesets = cubit.get_nodeset_count()
num_sidesets = cubit.get_sideset_count()

f = open("points.csv", "w")

##########################
########Header##############
##########################
#f.write("Header \n")
#f.write("%d %d %d \n", % (num_nodes, num_elements, num_edges) )

##########################
#####Write Points##############
##########################

f.write("Coordinates\n")
f.write("%d\n" % (num_nodes))
for i in range(1, num_nodes+1):
   xx = cubit.get_nodal_coordinates(i)
   x = xx[0]
   y = xx[1]
   z = xx[2]
#   print x
#   print y
   f.write( "%f %f %f\n" % (x, y, z) )


#########################
######Write Quads############
#########################
f.write("Elements\n")
f.write("%d\n"%  (num_elements))
for i in range(1, num_elements+1):
   nodes = cubit.get_connectivity('face',i)
   f.write("%d %d %d %d\n" % (nodes[0]-1, nodes[1]-1, nodes[2]-1, nodes[3]-1))

#########################
######Write Edges############
#########################
f.write("Edges\n")
f.write("%d\n" % num_edges)
for i in range(1, num_edges+1):
   nodes = cubit.get_connectivity('edge', i)
   f.write("%d %d\n" % (nodes[0]-1, nodes[1]-1))

f.write("Boundary Conditions\n")
#########################
#####Write Nodesets###########
#########################
f.write("Nodesets\n")
f.write("%d\n" % num_nodesets)
for i in range(1, num_nodesets+1):
   num_nodes_nodeset = cubit.get_nodeset_node_count(i)
   list = cubit.get_nodeset_nodes_inclusive(i)
   f.write("%d\n" % num_nodes_nodeset)
   for j in list:
      f.write("%d\n" % (j - 1))

#########################
#####Write Sidesets###########
#########################
f.write("Sidesets\n")
f.write("%d\n" % num_sidesets)
id_list = cubit.get_sideset_id_list()
for i in id_list:
   edge_list = cubit.get_sideset_edges(i)
   num_edges_sideset = len(edge_list)
   f.write("%d\n" % num_edges_sideset)
   for j in edge_list:
      f.write("%d\n" % (j-1))

f.close()


