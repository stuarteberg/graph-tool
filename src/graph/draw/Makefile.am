## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = $(MOD_CPPFLAGS) $(CAIROMM_CFLAGS)

AM_CFLAGS = $(AM_CXXFLAGS)

libgraph_tool_drawdir = $(MOD_DIR)/draw

libgraph_tool_draw_LTLIBRARIES = libgraph_tool_draw.la

libgraph_tool_draw_la_includedir = $(MOD_DIR)/include/draw

libgraph_tool_draw_la_LIBADD = $(MOD_LIBADD)

libgraph_tool_draw_la_LDFLAGS = $(MOD_LDFLAGS) $(CAIROMM_LIBS)

libgraph_tool_draw_la_SOURCES = \
    graph_cairo_draw.cc \
    graph_tree_cts.cc

libgraph_tool_draw_la_include_HEADERS = 

libgt_pycairo_auxdir = $(MOD_DIR)/draw

libgt_pycairo_aux_LTLIBRARIES = libgt_pycairo_aux.la

libgt_pycairo_aux_la_includedir = $(MOD_DIR)/include/draw

libgt_pycairo_aux_la_LIBADD = $(MOD_LIBADD)

libgt_pycairo_aux_la_LDFLAGS = $(MOD_LDFLAGS) $(CAIROMM_LIBS)

libgt_pycairo_aux_la_SOURCES = \
    pycairo_aux.cc

libgt_pycairo_aux_la_include_HEADERS = 
