#include <lua.h>
#include <lauxlib.h>
#include <lualib.h>
#include "micro-app-cuda.h"

int graph_init_aos(char* graph_type, int npoints, int nedges,
        struct edge* edges, char* fname) {
    lua_State *L;
    int i, j, k, v;

    L = luaL_newstate();
    luaL_openlibs(L);
    luaL_loadfile(L, "graph.lua");
    lua_pcall(L, 0, 0, 0);

    lua_getglobal(L, "create_graph");
    lua_pushstring(L, graph_type);
    lua_pushinteger(L, NPOINTS);
    lua_pushinteger(L, NEDGES);

    if (fname != NULL) {
        lua_pushstring(L, fname);
        lua_call(L, 4, 1);
    } else {
        lua_call(L, 3, 1);
    }

    /* Table is now sitting at the top of the stack */
    i = 0;
    lua_pushnil(L);  /* Make sure lua_next starts at beginning */
    while (lua_next(L, -2) != 0) {
        // fetch first key
        k = lua_tointeger(L, -2);
        lua_pushnil(L);
        while (lua_next(L, -2) != 0) { // loop over neighbors
            lua_pop(L,1);
            v = lua_tointeger(L, -1);

            // build edges array here
            edges[i].v0 = k - 1;
            edges[i].v1 = v - 1;

            for (j = 0; j < 3; j++) {
                edges[i].v0_pt_data[j] = 0;
                edges[i].v1_pt_data[j] = 0;
            }

            i++;
        }
        lua_pop(L,1);
    }

    lua_close(L);

    if (i == 0) {
        return -1;
    } else {
        return 0;
    }
}

int graph_init_soa(char* graph_type, int npoints, int nedges,
        struct graph* gr, char* fname) {
    lua_State *L;
    int i, j, k, v;

    L = luaL_newstate();
    luaL_openlibs(L);
    luaL_loadfile(L, "graph.lua");
    lua_pcall(L, 0, 0, 0);

    lua_getglobal(L, "create_graph");
    lua_pushstring(L, graph_type);
    lua_pushinteger(L, NPOINTS);
    lua_pushinteger(L, NEDGES);

    if (fname != NULL) {
        lua_pushstring(L, fname);
        lua_call(L, 4, 1);
    } else {
        lua_call(L, 3, 1);
    }

    /* Table is now sitting at the top of the stack */
    i = 0;
    lua_pushnil(L);  /* Make sure lua_next starts at beginning */
    while (lua_next(L, -2) != 0) {
        // fetch first key
        k = lua_tointeger(L, -2);
        lua_pushnil(L);
        while (lua_next(L, -2) != 0) { // loop over neighbors
            lua_pop(L,1);
            v = lua_tointeger(L, -1);

            // build the graph
            gr->v0[i] = k - 1;
            gr->v1[i] = v - 1;

            for (j = 0; j < 3; j++) {
                gr->v0_data[i][j] = 0;
                gr->v1_data[i][j] = 0;
            }

            i++;
        }
        lua_pop(L,1);
    }

    lua_close(L);

    if (i == 0) {
        return -1;
    } else {
        return 0;
    }
}
