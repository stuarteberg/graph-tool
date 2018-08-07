## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = $(MOD_CPPFLAGS)

AM_CFLAGS = $(AM_CXXFLAGS)

libgraph_tool_dynamicsdir = $(MOD_DIR)/dynamics

libgraph_tool_dynamics_LTLIBRARIES = libgraph_tool_dynamics.la

libgraph_tool_dynamics_la_includedir = $(MOD_DIR)/include/dynamics

libgraph_tool_dynamics_la_LIBADD = $(MOD_LIBADD)

libgraph_tool_dynamics_la_LDFLAGS = $(MOD_LDFLAGS)

libgraph_tool_dynamics_la_SOURCES = \
    graph_continuous.cc \
    graph_discrete.cc \
    graph_dynamics.cc

libgraph_tool_dynamics_la_include_HEADERS = \
    graph_continuous.hh \
    graph_discrete.hh

