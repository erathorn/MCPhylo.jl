##iterator function:
#input is input tree
#list of nodes used(LoNU) = []
#list of trees(LoT) = [input tree]
#while LoT not empty:
#curtree = LoT.pop
#copyoftree = curtree
#for node x in [all nodes in copyoftree]
#dest = random node from copyoftree (maybe keep track of this somehow 4speed)
#newtree = transformfunction(copyoftree,x, dest)
#copyoftree = curtree
# if newtree.length() < curtree.length()
#LoT.add(newtree)
#
#outside while loop: return curtree
#


##tree transforming function:
#input is root tree, thing to move(TtM), and destination(dest)
# removechild(mother of TtM, TtM)
#IF BINARY:
# make new node, addchild(newnode, TtM)
# addchild(newnode, dest)
#addchild(dest.mother, newnode)
# removechild(dest.mother,dest)
#IF NOT BINARY:
# addchild(dest.mother,TtM)
#

function transform(root::T, TtM::T, dest::T) where T<:AbstractNode
    motherofTtM = TtM.mother
    
