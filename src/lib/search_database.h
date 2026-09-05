// A store of already-visited points in conformation space.
//
// Stock Vina does not use one: bfgs() always runs the local optimisation. The
// QuickVina family's contribution is to consult such a store first and skip the
// local search when a first-order extrapolation from a nearby visited point
// already predicts what would be found. This header declares only the
// interface, so that Vina itself carries no policy and no cost.
//
// The interface is QuickVina-W's, unchanged:
//
//     int  interesting(conf x, double f, change g, int excluded)
//     void add(conf x, double f, change g)
//
// `interesting` returns a NEGATIVE value if the point is worth optimising, and
// a NON-NEGATIVE one if it is not -- in which case the value is also the index
// of the entry that matched, which the caller threads into the next store's
// `excluded` argument so that one entry cannot satisfy both checks.
//
// Passing a null pointer disables a store, and two null stores are exactly
// stock Vina. QuickVina-W already relied on that, guarding its consult with
// `if (tried)`, and already used one interface for both of its stores; this
// only gives them the common base that its own source anticipated, in the
// commented-out base class on `class circularvisited`.

#ifndef VINA_SEARCH_DATABASE_H
#define VINA_SEARCH_DATABASE_H

#include "conf.h"

struct search_database {
	virtual ~search_database() {}

	// Arguments are taken by value, matching QuickVina's own signature.
	virtual int interesting(conf x, double f, change g, int excluded) = 0;
	virtual void add(conf x, double f, change g) = 0;
};

#endif
