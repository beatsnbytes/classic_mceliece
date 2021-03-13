kat_kem.rsp: kat
	./run

kat: Makefile nist/kat_kem.c nist/rng.c nist/rng.h randombytes.h benes.c bm.c controlbits.c decrypt.c encrypt.c fft.c fft_tr.c gf.c operations.c pk_gen.c sk_gen.c vec.c 
ifeq ($(ARCH), arm)
	./build.aarch64 $(OPT) $(VEC) $(KAT)
else
	./build $(OPT) $(VEC) $(KAT)
endif

