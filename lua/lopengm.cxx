
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
typedef opengm::AStar<GraphicalModel, opengm::Minimizer> ASTAR;

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
  // args:
  int table_states = 1;
  int table_args = 2;
  int table_energies = 3;

  // get sizes
  size_t nb_variables = lua_objlen(L, 1);
  size_t nb_factors = lua_objlen(L, 2);
  size_t nb_energies = lua_objlen(L, 2);
  assert(nb_energies == nb_factors);

  // create state space
  std::vector<size_t> numbersOfStates(nb_variables, 2);
  for (size_t i=0; i<nb_variables; ++i) {
    lua_rawgeti(L, table_states, i+1);
    numbersOfStates[i] = lua_tonumber(L, -1);
    lua_pop(L, 1);
  }
  Space *space = new Space(numbersOfStates.begin(), numbersOfStates.end());

  // create graphical model
  GraphicalModel *gm = new GraphicalModel();

  // register all factors into graphical model
  for (size_t i=0; i<nb_factors; ++i) {
    // args (2 max for now)
    size_t args[2];

    // get next arg list:
    lua_rawgeti(L, table_args, i+1); size_t nargs = lua_objlen(L, -1);
    assert(nargs<=2);

    // get argument index
    for (size_t l=0; l<nargs; ++l) {
      lua_rawgeti(L, -1, l+1); args[l] = lua_tonumber(L,-1) - 1; lua_pop(L,1);
    }
    lua_pop(L, 1);

    // declare factor
    Factor factor(*space, args, args+nargs);
  
    // get next factor:
    lua_rawgeti(L, table_energies, i+1);

    // set energy for each state
    int kkk = 1;
    if (nargs == 1) {
      for (size_t k=0; k<numbersOfStates[args[0]]; ++k) {
        lua_rawgeti(L, -1, kkk++);
        factor(k) = (Energy)lua_tonumber(L,-1); 
        lua_pop(L,1);
      } 
    } else if (nargs == 2) {
      for (size_t k=0; k<numbersOfStates[args[0]]; ++k) {
        for (size_t l=0; l<numbersOfStates[args[1]]; ++l) {
          lua_rawgeti(L, -1, kkk++);
          factor(k,l) = (Energy)lua_tonumber(L,-1); 
          lua_pop(L,1);
        }
      }
    }
    lua_pop(L, 1);

    // register factor
    gm->addFactor(factor);
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
  const char *modec = lua_tostring(L,2);
  int verbose = lua_toboolean(L,3);
  if (lua_isnumber(L,2)) max_steps = lua_tonumber(L,4);
  std::string mode(modec);

  // state is used to hold the result of optimization
  State *state = new State(gm->space().dimension());

  // algorithm select
  if (mode.compare("bp") == 0) {
    // setup Belief Propagation
    BP::Parameter para;
    para.maximumNumberOfSteps_ = max_steps;
    para.damping_ = 0;
    BP bp(*gm, para);
    if(verbose) {
      std::cout << "<opengm> optimizing using Belief Propagation " << std::endl;
      opengm::BeliefPropagationVerboseVisitor<BP> visitor;
      bp.infer(visitor);
    } else {
      bp.infer();
    }
    bp.arg(*state);

  } else if (mode.compare("a*") == 0) {
    // setup A*
    ASTAR::Parameter para;
    ASTAR astar(*gm, para);
    if(verbose) {
      std::cout << "<opengm> optimizing using A* " << std::endl;
      opengm::AStarVisitor<ASTAR,true> visitor;
      astar.infer(visitor);
    } else {
      astar.infer();
    }
    astar.arg(*state);

  } else {
    THError("<opengm.Graph.optimize> method must be one of: a* | bp");
  }

  // save optimal state
  if (g->state) delete g->state;
  g->state = state;

  // done
  lua_settop(L,1);
  return 1;
}

static int Graph_states (lua_State *L)
{
  // args
  Graph *g = lua_checkGraph(L, 1);

  if (g->state == NULL) {
    printf("<opengm> states are not available yet, call graph:optimze() first\n");
    return 0;
  }

  // fill tensor with current config
  THDoubleTensor *states = THDoubleTensor_newWithSize1d(g->space->dimension());
  double *data = THDoubleTensor_data(states);
  for(int y=0; y<g->state->size(); y++) {
    data[y] = g->state->at(y);
  }

  // return tensor
  luaT_pushudata(L, states, luaT_checktypename2id(L, "torch.DoubleTensor"));
  return 1;
}

static const luaL_reg Graph_methods[] = {
  {"new",      Graph_new},
  {"optimize", Graph_optimize},
  {"states",   Graph_states},
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
