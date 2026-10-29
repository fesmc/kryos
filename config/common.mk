# Shared build configuration for kryos (dependency wiring).
#
# Loaded *after* the compiler and machine fragments (configme assembles them in
# the order: compiler -> machine -> netCDF -> common). This file references
# variables those provide: FFLAGS / FFLAGS_OPENMP (compiler) and LIB_NC
# (machine or auto-detected netCDF). It defines variables only; the rules that
# use them are in config/Makefile.
#
# Every dependency is a sibling checkout (or symlink) at the repo root, as in
# yelmox: fesm-utils, yelmo, FastIsostasy, chion, insol, ...

# fesm-utils: the only mandatory dependency (serial build by default here;
# swapped to the OpenMP variant below).
FESMUTILSROOT = fesm-utils
FESMUTILSLIBDIR = ${FESMUTILSROOT}/include-serial

# fesm-utils' external libraries, needed on the link line by yelmo (LIS) and by
# FastHydrology/FastIsostasy (FFTW).
LISROOT  = ${FESMUTILSROOT}/lis/lis-serial
FFTWROOT = ${FESMUTILSROOT}/fftw/fftw-serial
LIB_FFTW_LIBS = -lfftw3 -lm

# OpenMP build (make openmp=1): swap the serial dependency builds for their
# OpenMP variants and append the compiler's OpenMP flag (FFLAGS_OPENMP, set in
# the compiler fragment). The host must link the variant matching its build.
ifeq ($(openmp), 1)
    FESMUTILSLIBDIR = ${FESMUTILSROOT}/include-omp
    LISROOT  = ${FESMUTILSROOT}/lis/lis-omp
    FFTWROOT = ${FESMUTILSROOT}/fftw/fftw-omp
    LIB_FFTW_LIBS = -lfftw3_omp -lfftw3 -lm

    FFLAGS += $(FFLAGS_OPENMP)
endif

# Position-independent code (make pic=1): needed to link the objects into a
# shared library (a future C API); the static deps must be built with pic=1 too.
ifeq ($(pic), 1)
    FFLAGS += -fPIC
endif

INC_FESMUTILS = -I${FESMUTILSLIBDIR}
LIB_FESMUTILS = -L${FESMUTILSLIBDIR} -lfesmutils
FESMUTILS_ARCHIVE = $(FESMUTILSLIBDIR)/libfesmutils.a

INC_LIS = -I${LISROOT}/include
LIB_LIS = -L${LISROOT}/lib -llis
LIS_ARCHIVE = $(LISROOT)/lib/liblis.a

INC_FFTW = -I${FFTWROOT}/include
LIB_FFTW = -L${FFTWROOT}/lib $(LIB_FFTW_LIBS)
FFTW_ARCHIVE = $(FFTWROOT)/lib/libfftw3.a

# Backend library roots (used only when the backend is switched on below).
YELMOROOT = yelmo
INC_YELMO = -I${YELMOROOT}/libyelmo/include
LIB_YELMO = -L${YELMOROOT}/libyelmo/include -lyelmo
YELMO_ARCHIVE = $(YELMOROOT)/libyelmo/include/libyelmo.a

# yelmo's own dependencies: its static archive does not pull them in
# transitively, so the link line names them (as in yelmox).
LIB_FASTHYDRO = -L${YELMOROOT}/FastHydrology/include -lfasthydro
LIB_ELSA      = -L${YELMOROOT}/elsa/libelsa/include -lelsa
LIB_TRACER    = -L${YELMOROOT}/tracer/libtracer/include -ltracer

ISOSTASYROOT = FastIsostasy
INC_ISOSTASY = -I${ISOSTASYROOT}/libisostasy/include
LIB_ISOSTASY = -L${ISOSTASYROOT}/libisostasy/include -lisostasy
ISOSTASY_ARCHIVE = $(ISOSTASYROOT)/libisostasy/include/libisostasy.a

CHIONROOT = chion
INC_CHION = -I${CHIONROOT}/libchion/include
LIB_CHION = -L${CHIONROOT}/libchion/include -lchion
CHION_ARCHIVE = $(CHIONROOT)/libchion/include/libchion.a

# insol: serial, no dependencies.
INSOLROOT = insol
INC_INSOL = -I${INSOLROOT}/libinsol/include
LIB_INSOL = -L${INSOLROOT}/libinsol/include -linsol
INSOL_ARCHIVE = $(INSOLROOT)/libinsol/include/libinsol.a

# The roots of the input/ copies (input/sync.txt), for make sync-input and the
# staleness check of make check.
input_roots = FESMUTILSROOT=$(FESMUTILSROOT) YELMOROOT=$(YELMOROOT) CHIONROOT=$(CHIONROOT) \
              INSOLROOT=$(INSOLROOT) VILMAROOT=$(VILMAROOT)

REMBOROOT = rembo1
INC_REMBO = -I${REMBOROOT}/librembo/include
LIB_REMBO = -L${REMBOROOT}/librembo/include -lrembo
REMBO_ARCHIVE = $(REMBOROOT)/librembo/include/librembo.a

# VILMA2 (geo = vilma2): libvilma, and the SHTns that fesm-utils builds (make
# shtns, once; the serial or OpenMP variant as for LIS and FFTW). VILMA2 must be
# built with the same openmp= as kryos and fesm-utils.
VILMAROOT = vilma2
INC_VILMA = -I${VILMAROOT}/libvilma/include
LIB_VILMA = -L${VILMAROOT}/libvilma/include -lvilma
VILMA_ARCHIVE = $(VILMAROOT)/libvilma/include/libvilma.a

SHTNS_NAME = shtns
SHTNS_DIR  = shtns-serial
ifeq ($(openmp), 1)
    SHTNS_NAME = shtns_omp
    SHTNS_DIR  = shtns-omp
endif
LIB_SHTNS     = -L${FESMUTILSROOT}/SHTns/$(SHTNS_DIR)/lib -l$(SHTNS_NAME)
SHTNS_ARCHIVE = ${FESMUTILSROOT}/SHTns/$(SHTNS_DIR)/lib/lib$(SHTNS_NAME).a

# Backend switches (design section 9). Every backend library is optional: on,
# its INC/LIB/archive are added and its adapters are compiled; off, each
# adapter module is provided by its same-interface stub
# (src/backends/<component>/kryos_<component>_<backend>_stub.f90), whose init
# stops with "rebuild with <backend>=1". The adapter is selected by setting
# <module>_obj, which config/deps.mk requires for every module defined twice.
yelmo        ?= 1
fastisostasy ?= 1
chion        ?= 1
insol        ?= 1
cism         ?= 0
pism         ?= 0
vilma        ?= 0
rembo        ?= 0

# The switch settings, recorded in the object-dir stamp (config/Makefile).
backend_config = yelmo$(yelmo)-fastisostasy$(fastisostasy)-chion$(chion)-insol$(insol)-cism$(cism)-pism$(pism)-vilma$(vilma)-rembo$(rembo)

INC_BACKENDS   =
LIB_BACKENDS   =
dep_archives   = $(FESMUTILS_ARCHIVE)
backend_roots  =
backend_static =
test_skip      =
need_lis_fftw  = 0

# The libraries whose version the build records (kryos_build_info, for the
# restart manifest): <name>:<root>, named as the build switch.
build_info_libs = fesm-utils:$(FESMUTILSROOT)

ifeq ($(yelmo), 1)
    kryos_ice_yelmo_obj = $(objdir)/kryos_ice_yelmo.o
    kryos_exp_relax_obj = $(objdir)/kryos_exp_relax.o
    kryos_exp_opt_cb_ref_obj = $(objdir)/kryos_exp_opt_cb_ref.o
    kryos_exp_opt_tf_obj = $(objdir)/kryos_exp_opt_tf.o
    kryos_exp_recon_obj = $(objdir)/kryos_exp_recon.o
    INC_BACKENDS   += $(INC_YELMO)
    LIB_BACKENDS   += $(LIB_YELMO) $(LIB_FASTHYDRO) $(LIB_ELSA) $(LIB_TRACER)
    dep_archives   += $(YELMO_ARCHIVE)
    backend_roots  += $(YELMOROOT)
    build_info_libs += yelmo:$(YELMOROOT)
    backend_static += yelmo-static
    need_lis_fftw   = 1
else
    kryos_ice_yelmo_obj = $(objdir)/kryos_ice_yelmo_stub.o
    kryos_exp_relax_obj = $(objdir)/kryos_exp_relax_stub.o
    kryos_exp_opt_cb_ref_obj = $(objdir)/kryos_exp_opt_cb_ref_stub.o
    kryos_exp_opt_tf_obj = $(objdir)/kryos_exp_opt_tf_stub.o
    kryos_exp_recon_obj = $(objdir)/kryos_exp_recon_stub.o
    test_skip += $(testdir)/test_run_sequence.f90 $(testdir)/test_exp_opt_tf.f90
endif

ifeq ($(fastisostasy), 1)
    kryos_geo_fastisostasy_obj = $(objdir)/kryos_geo_fastisostasy.o
    INC_BACKENDS   += $(INC_ISOSTASY)
    LIB_BACKENDS   += $(LIB_ISOSTASY)
    dep_archives   += $(ISOSTASY_ARCHIVE)
    backend_roots  += $(ISOSTASYROOT)
    build_info_libs += fastisostasy:$(ISOSTASYROOT)
    backend_static += isostasy-static
    need_lis_fftw   = 1
else
    kryos_geo_fastisostasy_obj = $(objdir)/kryos_geo_fastisostasy_stub.o
    test_skip += $(testdir)/test_geo_fastisostasy.f90 $(testdir)/test_globe.f90 $(testdir)/test_run_sequence.f90
endif

ifeq ($(chion), 1)
    # The chion host (kryos_surface_chion.f90) takes its daily insolation from insol.
    ifneq ($(insol), 1)
        $(error chion=1 needs insol=1: the chion host uses libinsol)
    endif
    kryos_surface_chion_obj = $(objdir)/kryos_surface_chion.o
    INC_BACKENDS   += $(INC_CHION)
    LIB_BACKENDS   += $(LIB_CHION)
    dep_archives   += $(CHION_ARCHIVE)
    backend_roots  += $(CHIONROOT)
    build_info_libs += chion:$(CHIONROOT)
    backend_static += chion-static
else
    kryos_surface_chion_obj = $(objdir)/kryos_surface_chion_stub.o
    test_skip += $(testdir)/test_run_sequence.f90
endif

ifeq ($(insol), 1)
    INC_BACKENDS   += $(INC_INSOL)
    LIB_BACKENDS   += $(LIB_INSOL)
    dep_archives   += $(INSOL_ARCHIVE)
    backend_roots  += $(INSOLROOT)
    build_info_libs += insol:$(INSOLROOT)
    backend_static += insol-static
endif

ifeq ($(vilma), 1)
    kryos_geo_vilma2_obj = $(objdir)/kryos_geo_vilma2.o
    INC_BACKENDS   += $(INC_VILMA)
    LIB_BACKENDS   += $(LIB_VILMA)
    dep_archives   += $(VILMA_ARCHIVE) $(SHTNS_ARCHIVE)
    backend_roots  += $(VILMAROOT)
    build_info_libs += vilma:$(VILMAROOT)
    backend_static += vilma-static
    need_lis_fftw   = 1
    LIB_SHTNS_USED  = $(LIB_SHTNS)
    test_skip += $(testdir)/test_geo_vilma2_stub.f90
else
    kryos_geo_vilma2_obj = $(objdir)/kryos_geo_vilma2_stub.o
    test_skip += $(testdir)/test_geo_vilma2.f90
endif

# Backends whose adapters arrive in later phases (design section 14).
ifeq ($(rembo), 1)
    $(error rembo=1: the climate/rembo adapter does not exist yet (phase 2))
endif
ifeq ($(cism), 1)
    $(error cism=1: the ice/cism adapter does not exist yet (phase 5))
endif
ifeq ($(pism), 1)
    $(error pism=1: the ice/pism adapter does not exist yet (phase 5))
endif

# fesm-utils' LIS and FFTW, when a backend needs them.
LIB_EXTERNAL =
LIB_SHTNS_USED ?=
ifeq ($(need_lis_fftw), 1)
    LIB_EXTERNAL  = $(LIB_LIS) $(LIB_SHTNS_USED) $(LIB_FFTW)
    dep_archives += $(LIS_ARCHIVE) $(FFTW_ARCHIVE)
endif

# Phase 3 WP1: the two-domain test needs what the sequence test needs.
ifneq ($(findstring test_run_sequence.f90,$(test_skip)),)
    test_skip += $(testdir)/test_run_domains.f90
endif

# Extra link flags. -Wl,-zmuldefs works around duplicate symbols in the static
# deps (the default on Linux). A machine fragment can disable it by setting
# `LFLAGS_EXTRA =` (macOS ld rejects -zmuldefs, so the macbook fragment does).
LFLAGS_EXTRA ?= -Wl,-zmuldefs

# Backend libraries before fesm-utils, which they all use; netCDF and the
# external libraries last.
LFLAGS = $(LIB_BACKENDS) $(LIB_FESMUTILS) $(LIB_NC) $(LIB_EXTERNAL) $(LFLAGS_EXTRA)
