include ../../makeconf

INCPATH = -I../../include
CFLAGS = $(MOREFLAGS) $(WARNLEVEL) -pedantic $(INCPATH) -I../RANGECODER

SOURCES.C    =  point.cpp\
		vertex.cpp\
		edge.cpp\
		triangle.cpp\
		tin.cpp\
		io.cpp\
		checkAndRepair.cpp
		
OBJECTS.O    =  $(SOURCES.C:.cpp=.o)

../OBJECTS/mesh_structure.o: $(OBJECTS.O)
	$(LD) -r -o ../OBJECTS/mesh_structure.o $(OBJECTS.O)

depend:
	makedepend -- $(CFLAGS) -- $(SOURCES.C)
	$(RM) -f Makefile.bak

clean :
	$(RM) -f *.o ../OBJECTS/mesh_structure.o

.cpp.o:
	$(CXX) $(CFLAGS) -c -o $@ $<

# DO NOT DELETE

point.o: ../../include/point.h ../../include/j_mesh.h /usr/include/stdio.h
point.o: /usr/include/features.h /usr/include/sys/cdefs.h
point.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
point.o: /usr/include/gnu/stubs-32.h /usr/include/bits/types.h
point.o: /usr/include/bits/typesizes.h /usr/include/libio.h
point.o: /usr/include/_G_config.h /usr/include/wchar.h
point.o: /usr/include/bits/wchar.h /usr/include/gconv.h
point.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
point.o: /usr/include/math.h /usr/include/bits/huge_val.h
point.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
point.o: /usr/include/limits.h /usr/include/bits/posix1_lim.h
point.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
point.o: /usr/include/bits/posix2_lim.h /usr/include/stdlib.h
point.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/endian.h
point.o: /usr/include/bits/endian.h /usr/include/sys/select.h
point.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
point.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
point.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
point.o: /usr/include/errno.h /usr/include/bits/errno.h
point.o: /usr/include/linux/errno.h /usr/include/asm/errno.h
point.o: /usr/include/asm-generic/errno.h
point.o: /usr/include/asm-generic/errno-base.h
vertex.o: ../../include/vertex.h ../../include/j_mesh.h /usr/include/stdio.h
vertex.o: /usr/include/features.h /usr/include/sys/cdefs.h
vertex.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
vertex.o: /usr/include/gnu/stubs-32.h /usr/include/bits/types.h
vertex.o: /usr/include/bits/typesizes.h /usr/include/libio.h
vertex.o: /usr/include/_G_config.h /usr/include/wchar.h
vertex.o: /usr/include/bits/wchar.h /usr/include/gconv.h
vertex.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
vertex.o: /usr/include/math.h /usr/include/bits/huge_val.h
vertex.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
vertex.o: /usr/include/limits.h /usr/include/bits/posix1_lim.h
vertex.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
vertex.o: /usr/include/bits/posix2_lim.h ../../include/list.h
vertex.o: ../../include/point.h ../../include/edge.h ../../include/triangle.h
vertex.o: ../../include/jmesh.h ../../include/matrix.h ../../include/tin.h
vertex.o: /usr/include/stdlib.h /usr/include/sys/types.h /usr/include/time.h
vertex.o: /usr/include/endian.h /usr/include/bits/endian.h
vertex.o: /usr/include/sys/select.h /usr/include/bits/select.h
vertex.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
vertex.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
vertex.o: /usr/include/alloca.h /usr/include/errno.h
vertex.o: /usr/include/bits/errno.h /usr/include/linux/errno.h
vertex.o: /usr/include/asm/errno.h /usr/include/asm-generic/errno.h
vertex.o: /usr/include/asm-generic/errno-base.h
edge.o: ../../include/edge.h ../../include/j_mesh.h /usr/include/stdio.h
edge.o: /usr/include/features.h /usr/include/sys/cdefs.h
edge.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
edge.o: /usr/include/gnu/stubs-32.h /usr/include/bits/types.h
edge.o: /usr/include/bits/typesizes.h /usr/include/libio.h
edge.o: /usr/include/_G_config.h /usr/include/wchar.h
edge.o: /usr/include/bits/wchar.h /usr/include/gconv.h
edge.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
edge.o: /usr/include/math.h /usr/include/bits/huge_val.h
edge.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
edge.o: /usr/include/limits.h /usr/include/bits/posix1_lim.h
edge.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
edge.o: /usr/include/bits/posix2_lim.h ../../include/vertex.h
edge.o: ../../include/list.h ../../include/point.h ../../include/triangle.h
edge.o: ../../include/jmesh.h ../../include/matrix.h ../../include/tin.h
triangle.o: ../../include/triangle.h ../../include/jmesh.h
triangle.o: ../../include/j_mesh.h /usr/include/stdio.h
triangle.o: /usr/include/features.h /usr/include/sys/cdefs.h
triangle.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
triangle.o: /usr/include/gnu/stubs-32.h /usr/include/bits/types.h
triangle.o: /usr/include/bits/typesizes.h /usr/include/libio.h
triangle.o: /usr/include/_G_config.h /usr/include/wchar.h
triangle.o: /usr/include/bits/wchar.h /usr/include/gconv.h
triangle.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
triangle.o: /usr/include/math.h /usr/include/bits/huge_val.h
triangle.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
triangle.o: /usr/include/limits.h /usr/include/bits/posix1_lim.h
triangle.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
triangle.o: /usr/include/bits/posix2_lim.h ../../include/list.h
triangle.o: ../../include/point.h ../../include/vertex.h ../../include/edge.h
triangle.o: ../../include/matrix.h ../../include/tin.h /usr/include/stdlib.h
triangle.o: /usr/include/sys/types.h /usr/include/time.h
triangle.o: /usr/include/endian.h /usr/include/bits/endian.h
triangle.o: /usr/include/sys/select.h /usr/include/bits/select.h
triangle.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
triangle.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
triangle.o: /usr/include/alloca.h
tin.o: ../../include/tin.h ../../include/jmesh.h ../../include/j_mesh.h
tin.o: /usr/include/stdio.h /usr/include/features.h /usr/include/sys/cdefs.h
tin.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
tin.o: /usr/include/gnu/stubs-32.h /usr/include/bits/types.h
tin.o: /usr/include/bits/typesizes.h /usr/include/libio.h
tin.o: /usr/include/_G_config.h /usr/include/wchar.h
tin.o: /usr/include/bits/wchar.h /usr/include/gconv.h
tin.o: /usr/include/bits/stdio_lim.h /usr/include/bits/sys_errlist.h
tin.o: /usr/include/math.h /usr/include/bits/huge_val.h
tin.o: /usr/include/bits/mathdef.h /usr/include/bits/mathcalls.h
tin.o: /usr/include/limits.h /usr/include/bits/posix1_lim.h
tin.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
tin.o: /usr/include/bits/posix2_lim.h ../../include/list.h
tin.o: ../../include/point.h ../../include/vertex.h ../../include/edge.h
tin.o: ../../include/triangle.h ../../include/matrix.h /usr/include/stdlib.h
tin.o: /usr/include/sys/types.h /usr/include/time.h /usr/include/endian.h
tin.o: /usr/include/bits/endian.h /usr/include/sys/select.h
tin.o: /usr/include/bits/select.h /usr/include/bits/sigset.h
tin.o: /usr/include/bits/time.h /usr/include/sys/sysmacros.h
tin.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
tin.o: /usr/include/string.h
io.o: ../../include/jmesh.h ../../include/j_mesh.h /usr/include/stdio.h
io.o: /usr/include/features.h /usr/include/sys/cdefs.h
io.o: /usr/include/bits/wordsize.h /usr/include/gnu/stubs.h
io.o: /usr/include/gnu/stubs-32.h /usr/include/bits/types.h
io.o: /usr/include/bits/typesizes.h /usr/include/libio.h
io.o: /usr/include/_G_config.h /usr/include/wchar.h /usr/include/bits/wchar.h
io.o: /usr/include/gconv.h /usr/include/bits/stdio_lim.h
io.o: /usr/include/bits/sys_errlist.h /usr/include/math.h
io.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
io.o: /usr/include/bits/mathcalls.h /usr/include/limits.h
io.o: /usr/include/bits/posix1_lim.h /usr/include/bits/local_lim.h
io.o: /usr/include/linux/limits.h /usr/include/bits/posix2_lim.h
io.o: ../../include/list.h ../../include/point.h ../../include/vertex.h
io.o: ../../include/edge.h ../../include/triangle.h ../../include/matrix.h
io.o: ../../include/tin.h /usr/include/stdlib.h /usr/include/sys/types.h
io.o: /usr/include/time.h /usr/include/endian.h /usr/include/bits/endian.h
io.o: /usr/include/sys/select.h /usr/include/bits/select.h
io.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
io.o: /usr/include/sys/sysmacros.h /usr/include/bits/pthreadtypes.h
io.o: /usr/include/alloca.h /usr/include/string.h /usr/include/ctype.h
checkAndRepair.o: ../../include/jmesh.h ../../include/j_mesh.h
checkAndRepair.o: /usr/include/stdio.h /usr/include/features.h
checkAndRepair.o: /usr/include/sys/cdefs.h /usr/include/bits/wordsize.h
checkAndRepair.o: /usr/include/gnu/stubs.h /usr/include/gnu/stubs-32.h
checkAndRepair.o: /usr/include/bits/types.h /usr/include/bits/typesizes.h
checkAndRepair.o: /usr/include/libio.h /usr/include/_G_config.h
checkAndRepair.o: /usr/include/wchar.h /usr/include/bits/wchar.h
checkAndRepair.o: /usr/include/gconv.h /usr/include/bits/stdio_lim.h
checkAndRepair.o: /usr/include/bits/sys_errlist.h /usr/include/math.h
checkAndRepair.o: /usr/include/bits/huge_val.h /usr/include/bits/mathdef.h
checkAndRepair.o: /usr/include/bits/mathcalls.h /usr/include/limits.h
checkAndRepair.o: /usr/include/bits/posix1_lim.h
checkAndRepair.o: /usr/include/bits/local_lim.h /usr/include/linux/limits.h
checkAndRepair.o: /usr/include/bits/posix2_lim.h ../../include/list.h
checkAndRepair.o: ../../include/point.h ../../include/vertex.h
checkAndRepair.o: ../../include/edge.h ../../include/triangle.h
checkAndRepair.o: ../../include/matrix.h ../../include/tin.h
checkAndRepair.o: ../../include/jqsort.h /usr/include/stdlib.h
checkAndRepair.o: /usr/include/sys/types.h /usr/include/time.h
checkAndRepair.o: /usr/include/endian.h /usr/include/bits/endian.h
checkAndRepair.o: /usr/include/sys/select.h /usr/include/bits/select.h
checkAndRepair.o: /usr/include/bits/sigset.h /usr/include/bits/time.h
checkAndRepair.o: /usr/include/sys/sysmacros.h
checkAndRepair.o: /usr/include/bits/pthreadtypes.h /usr/include/alloca.h
checkAndRepair.o: /usr/include/string.h
