CC       ?= cc
AR       ?= ar
CFLAGS   ?= -O3 -std=c99 -Wall -Wextra -Wpedantic
CPPFLAGS ?= -Iinclude
LDLIBS   ?= -lm
PREFIX   ?= /usr/local
INCDIR   ?= $(PREFIX)/include
LIBDIR   ?= $(PREFIX)/lib

BUILD = build
LIB_A = $(BUILD)/libamd_lite.a
LIB_SO = $(BUILD)/libamd_lite.so

all: $(LIB_A) $(LIB_SO)

$(BUILD):
	mkdir -p $(BUILD)

$(BUILD)/amd_lite.o: src/amd_lite.c src/amd_lite_impl.h include/amd_lite.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) -c src/amd_lite.c -o $@

$(BUILD)/amd_lite.pic.o: src/amd_lite.c src/amd_lite_impl.h include/amd_lite.h | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) -fPIC -c src/amd_lite.c -o $@

$(LIB_A): $(BUILD)/amd_lite.o
	$(AR) rcs $@ $<

$(LIB_SO): $(BUILD)/amd_lite.pic.o
	$(CC) -shared -o $@ $< $(LDLIBS)

$(BUILD)/test_amd_lite: tests/test_amd_lite.c $(LIB_A) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) tests/test_amd_lite.c $(LIB_A) $(LDLIBS) -o $@

$(BUILD)/amd_lite_demo: demo/amd_lite_demo.c $(LIB_A) | $(BUILD)
	$(CC) $(CPPFLAGS) $(CFLAGS) demo/amd_lite_demo.c $(LIB_A) $(LDLIBS) -o $@

test: $(BUILD)/test_amd_lite
	./$(BUILD)/test_amd_lite

demo: $(BUILD)/amd_lite_demo
	./$(BUILD)/amd_lite_demo 20

install: all
	mkdir -p $(DESTDIR)$(INCDIR) $(DESTDIR)$(LIBDIR)
	cp include/amd_lite.h $(DESTDIR)$(INCDIR)/
	cp $(LIB_A) $(LIB_SO) $(DESTDIR)$(LIBDIR)/

clean:
	rm -rf $(BUILD)

.PHONY: all test demo install clean
