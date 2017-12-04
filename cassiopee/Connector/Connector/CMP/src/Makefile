include ../Make.inc

OBJS = CMP_C.o  CMP.o  CMP_F.o CMP_For.o  CMP_FUtil.o

default: libBuffer.a


libBuffer.a: $(OBJS)
	ar rcs $@ $^

cleanall:
	@rm -f *.o *~ *.a *.mod
