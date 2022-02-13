include ../../makeconf

INCPATH = -I../../include -I../THIRD_PARTY
CFLAGS = $(MOREFLAGS) $(WARNLEVEL) -pedantic $(INCPATH)

SOURCES.C    =  binTree.cpp\
		dijkstraGraph.cpp\
		heap.cpp\
		list.cpp\
		clusterGraph.cpp\
		graph.cpp\
		jqsort.cpp\
		matrix.cpp

OBJECTS.O    =  $(SOURCES.C:.cpp=.o)

../OBJECTS/primitives.o: $(OBJECTS.O)
	$(LD) -r -o ../OBJECTS/primitives.o $(OBJECTS.O)

depend:
	makedepend -- $(CFLAGS) -- $(SOURCES.C)
	$(RM) -f Makefile.bak

clean :
	$(RM) -f *.o ../OBJECTS/primitives.o

.cpp.o:
	$(CXX) $(CFLAGS) -c -o $@ $<

# DO NOT DELETE

binTree.o: /usr/include/stdio.h /usr/include/features.h
binTree.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
binTree.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-32.h
binTree.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
binTree.o: /usr/include/libio.h /usr/include/_G_config.h /usr/include/wchar.h
binTree.o: /usr/include/bits/wchar.h /usr/include/gconv.h
binTree.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
binTree.o: ../../include/binTree.h
dijkstraGraph.o: ../../include/dijkstraGraph.h ../../include/heap.h
dijkstraGraph.o: ../../include/j_mesh.h /usr/include/stdio.h
dijkstraGraph.o: /usr/include/features.h /usr/include/sys/cdefs.h
dijkstraGraph.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
dijkstraGraph.o: /usr/include/gnu/stubs-32.h /usr/include/bits/types.h
dijkstraGraph.o: /usr/include/bits/typesizes.h /usr/include/libio.h
dijkstraGraph.o: /usr/include/_G_config.h /usr/include/wchar.h
dijkstraGraph.o: /usr/include/bits/wchar.h /usr/include/gconv.h
dijkstraGraph.o: /usr/include/bits/stdio_lim.h
dijkstraGraph.o: /usr/include/bits/sys_errlist.h /usr/include/math.h
dijkstraGraph.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
dijkstraGraph.o: /usr/include/bits/mathcalls.h /usr/include/limits.h
dijkstraGraph.o: /usr/include/bits/posix1_lim.h /usr/include/bits/local_lim.h
dijkstraGraph.o: /usr/include/linux/limits.h /usr/include/bits/posix2_lim.h
dijkstraGraph.o: ../../include/graph.h ../../include/list.h
heap.o: /usr/include/stdio.h /usr/include/features.h /usr/include/sys/cdefs.h
heap.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
heap.o: /usr/include/gnu/stubs-32.h /usr/include/bits/types.h
heap.o: /usr/include/bits/typesizes.h /usr/include/libio.h
heap.o: /usr/include/_G_config.h /usr/include/wchar.h
heap.o: /usr/include/bits/wchar.h /usr/include/gconv.h
heap.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
heap.o: /usr/include/stdlib.h /usr/include/sys/types.h /usr/include/time.h
heap.o: /usr/include/endian.h /usr/include/bits/endian.h
heap.o: /usr/include/sys/select.h /usr/include/bits/select.h
heap.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
heap.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
heap.o: /usr/include/alloca.h ../../include/heap.h ../../include/j_mesh.h
heap.o: /usr/include/math.h /usr/include/bits/huge_val.h
heap.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
heap.o: /usr/include/limits.h /usr/include/bits/posix1_lim.h
heap.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
heap.o: /usr/include/bits/posix2_lim.h
list.o: /usr/include/stdio.h /usr/include/features.h /usr/include/sys/cdefs.h
list.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
list.o: /usr/include/gnu/stubs-32.h /usr/include/bits/types.h
list.o: /usr/include/bits/typesizes.h /usr/include/libio.h
list.o: /usr/include/_G_config.h /usr/include/wchar.h
list.o: /usr/include/bits/wchar.h /usr/include/gconv.h
list.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
list.o: /usr/include/stdlib.h /usr/include/sys/types.h /usr/include/time.h
list.o: /usr/include/endian.h /usr/include/bits/endian.h
list.o: /usr/include/sys/select.h /usr/include/bits/select.h
list.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
list.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
list.o: /usr/include/alloca.h ../../include/list.h ../../include/jqsort.h
clusterGraph.o: ../../include/clusterGraph.h ../../include/heap.h
clusterGraph.o: ../../include/j_mesh.h /usr/include/stdio.h
clusterGraph.o: /usr/include/features.h /usr/include/sys/cdefs.h
clusterGraph.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
clusterGraph.o: /usr/include/gnu/stubs-32.h /usr/include/bits/types.h
clusterGraph.o: /usr/include/bits/typesizes.h /usr/include/libio.h
clusterGraph.o: /usr/include/_G_config.h /usr/include/wchar.h
clusterGraph.o: /usr/include/bits/wchar.h /usr/include/gconv.h
clusterGraph.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
clusterGraph.o: /usr/include/math.h /usr/include/bits/huge_val.h
clusterGraph.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
clusterGraph.o: /usr/include/limits.h /usr/include/bits/posix1_lim.h
clusterGraph.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
clusterGraph.o: /usr/include/bits/posix2_lim.h ../../include/graph.h
clusterGraph.o: ../../include/list.h
graph.o: ../../include/graph.h ../../include/list.h /usr/include/stdio.h
graph.o: /usr/include/features.h /usr/include/sys/cdefs.h
graph.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
graph.o: /usr/include/gnu/stubs-32.h /usr/include/bits/types.h
graph.o: /usr/include/bits/typesizes.h /usr/include/libio.h
graph.o: /usr/include/_G_config.h /usr/include/wchar.h
graph.o: /usr/include/bits/wchar.h /usr/include/gconv.h
graph.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
matrix.o: /usr/include/math.h /usr/include/features.h
matrix.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
matrix.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-32.h
matrix.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
matrix.o: /usr/include/bits/mathcalls.h ../../include/matrix.h
matrix.o: /usr/include/stdio.h /usr/include/bits/types.h
matrix.o: /usr/include/bits/typesizes.h /usr/include/libio.h
matrix.o: /usr/include/_G_config.h /usr/include/wchar.h
matrix.o: /usr/include/bits/wchar.h /usr/include/gconv.h
matrix.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
matrix.o: ../../include/list.h
