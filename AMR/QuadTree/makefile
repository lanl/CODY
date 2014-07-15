Main: quadTreeVis.o treeRenderer.o Neighbor.o OneLevel.o QuadTree.o
	g++ -pg -o vis -framework OpenGL -framework GLUT quadTreeVis.o QuadTree.o treeRenderer.o OneLevel.o Neighbor.o
quadTreeVis: quadTreeVis.cpp
	g++  -c -pg quadTreeVis.cpp -framework OpenGL -framework GLUT
treeRenderer: treeRenderer.cpp
	g++ -c -pg  treeRenderer.cpp -framework OpenGL -framework GLUT
Neighbor: Neighbor.cpp
	g++ -c -pg Neighbor.cpp
OneLevel: OneLevel.cpp
	g++ -c -pg OneLevel.cpp
QuadTree: QuadTree.cpp Application.h
	g++ -c -pg QuadTree.cpp

clean:
	rm vis quadTreeVis.o treeRenderer.o Neighbor.o OneLevel.o QuadTree.o

