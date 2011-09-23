
// standard libs
#include <string>
#include <sstream>

// opengm headers
#include <opengm/graphicalmodel.hxx>
#include <opengm/serialization.hxx>
#include <opengm/adder.hxx>
#include <opengm/minimizer.hxx>
#include <opengm/maxdistance.hxx>
#include <opengm/decomposer.hxx>
#include <opengm/inference/beliefpropagation.hxx>
#include <opengm/inference/treereweightedbeliefpropagation.hxx>
#include <opengm/inference/icm.hxx>
#include <opengm/inference/lazyflipper.hxx>
#include <opengm/inference/astar.hxx>

// Torch7 headers
#include <TH.h>
#include <luaT.h>

// define a couple of types that make life easier
typedef double Energy;
typedef opengm::DiscreteSpace Space;
typedef opengm::ExplicitFactor<Energy> Factor;
typedef opengm::GraphicalModel<Factor, opengm::Adder> GraphicalModel;
typedef std::vector<size_t> State;
typedef opengm::BeliefPropagation<GraphicalModel, opengm::Minimizer, opengm::MaxDistance> BP;

// class name
#define GRAPH "opengm.Graph"

// create a proper Lua class, to represent our graph + sample space
typedef struct {
  Space *space;
  GraphicalModel *model;
  State *state;
} Graph;

static Graph *lua_toGraph (lua_State *L, int index)
{
  Graph *graph = (Graph *)lua_touserdata(L, index);
  if (graph == NULL) luaL_typerror(L, index, GRAPH);
  return graph;
}

static Graph *lua_checkGraph (lua_State *L, int index)
{
  Graph *graph;
  luaL_checktype(L, index, LUA_TUSERDATA);
  graph = (Graph *)luaL_checkudata(L, index, GRAPH);
  if (graph == NULL) luaL_typerror(L, index, GRAPH);
  return graph;
}

static Graph *lua_pushGraph (lua_State *L)
{
  Graph *graph = (Graph *)lua_newuserdata(L, sizeof(Graph));
  graph->space = NULL;
  graph->model = NULL;
  graph->state = NULL;
  luaL_getmetatable(L, GRAPH);
  lua_setmetatable(L, -2);
  return graph;
}

static int Graph_new (lua_State *L)
{
  // args


  // discrete state space and graphical model
  size_t gridSize = 10;
  std::vector<size_t> numbersOfStates(gridSize*gridSize, 3);
  Space *space = new Space(numbersOfStates.begin(), numbersOfStates.end());
  GraphicalModel *gm = new GraphicalModel();

  // single site factors
  for(size_t j=0; j<space->dimension(); ++j) {
    Factor factor(*space, &j, &j+1);
    factor(0) = static_cast<Energy>(rand());
    factor(1) = static_cast<Energy>(RAND_MAX) - factor(0);
    factor(2) = static_cast<Energy>(rand());
    gm->addFactor(factor);
  }

  // Potts factors
  Energy alpha = 10000;
  for(size_t y0=0; y0<gridSize; ++y0) {
    for(size_t x0=0; x0<gridSize; ++x0) {
      if(x0 != gridSize-1) {
        size_t x1 = x0+1;
        size_t y1 = y0;
        size_t vi[] = {x0 + gridSize*y0, x1 + gridSize*y1};
        Factor factor(*space, vi, vi+2);
        factor(0,0) = 0.0;
        factor(0,1) = alpha;
        factor(1,0) = alpha;
        factor(1,1) = 0.0;
        factor(0,2) = 0.0;
        factor(2,0) = alpha;
        factor(1,2) = alpha;
        factor(2,1) = 0.0;
        gm->addFactor(factor);
      }
      if(y0 != gridSize-1) {
        size_t x1 = x0;
        size_t y1 = y0+1;
        size_t vi[] = {x0 + gridSize*y0, x1 + gridSize*y1};
        Factor factor(*space, vi, vi+2);
        factor(0,0) = 0.0;
        factor(0,1) = alpha;
        factor(1,0) = alpha;
        factor(1,1) = 0.0;
        factor(0,2) = 0.0;
        factor(2,0) = alpha;
        factor(1,2) = alpha;
        factor(2,1) = 0.0;
        gm->addFactor(factor);
      }
    }
  }

  // return graph
  Graph *graph = lua_pushGraph(L);
  graph->model = gm;
  graph->space = space;
  return 1;
}

static int Graph_optimize (lua_State *L)
{
  // args
  Graph *g = lua_checkGraph(L, 1);
  GraphicalModel *gm = g->model;
  int max_steps = 100;
  if (lua_isnumber(L,2)) max_steps = lua_tonumber(L,2);
  int verbose = lua_toboolean(L,3);

  // setup Belief Propagation
  State *state = new State(gm->space().dimension());
  BP::Parameter para;
  para.maximumNumberOfSteps_ = max_steps;
  para.damping_ = 0;
  BP bp(*gm, para);

  // optimize
  if(verbose) {
    std::cout << "<opengm> optimizing... " << std::endl;
    opengm::BeliefPropagationVerboseVisitor<BP> visitor;
    bp.infer(visitor);
  } else {
    bp.infer();
  }

  // save optimal state
  bp.arg(*state);
  if (g->state) delete g->state;
  g->state = state;

  // done
  return 0;
}

static const luaL_reg Graph_methods[] = {
  {"new",      Graph_new},
  {"optimize", Graph_optimize},
  {0, 0}
};

static int Graph_gc (lua_State *L)
{
  Graph *g = lua_toGraph(L, 1);
  if (g->space) delete g->space;
  if (g->model) delete g->model;
  if (g->state) delete g->state;
  return 0;
}

static int Graph_tostring (lua_State *L)
{
  Graph *g = lua_toGraph(L, 1);
  std::string str;
  str += "<opengm.Graph>";
  if (g->model) {
    std::stringstream s1;
    s1 << g->space->dimension();
    str += "\n  + nb of variables: " + s1.str();
    std::stringstream s2;
    s2 << g->model->numberOfFactors();
    str += "\n  + nb of factors: " + s2.str();
    if (g->model->isAcyclic())
      str += "\n  + graph is acyclic";
    else
      str += "\n  + graph has cycles";
  }
  if (g->state) {
    str += "\n  + current (optimized) variable states: ";
    str += "\n    ";
    int i = 0;
    for(int y=0; y<g->state->size(); y++) {
      if (((i%20)==0) && i>0) str += "\n    ";
      std::stringstream ss;
      ss << g->state->at(y);
      str += ss.str() + " ";
      i++;
    }
  }
  lua_pushfstring(L, "%s", str.c_str());
  return 1;
}

static const luaL_reg Graph_meta[] = {
  {"__gc",       Graph_gc},
  {"__tostring", Graph_tostring},
  {0, 0}
};

extern "C" {
  int luaopen_liblopengm (lua_State *L) {
    luaL_openlib(L, GRAPH, Graph_methods, 0);  /* create methods table,
                                                  add it to the globals */
    luaL_newmetatable(L, GRAPH);          /* create metatable for Graph,
                                             and add it to the Lua registry */
    luaL_openlib(L, 0, Graph_meta, 0);    /* fill metatable */
    lua_pushliteral(L, "__index");
    lua_pushvalue(L, -3);               /* dup methods table*/
    lua_rawset(L, -3);                  /* metatable.__index = methods */
    lua_pushliteral(L, "__metatable");
    lua_pushvalue(L, -3);               /* dup methods table*/
    lua_rawset(L, -3);                  /* hide metatable:
                                           metatable.__metatable = methods */
    lua_pop(L, 1);                      /* drop metatable */
    return 1;                           /* return methods on the stack */
  }
}
