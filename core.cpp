#include "includes.h"

class Simulation
{
	public:
		vector<ObjectSpec *> objectParts; // What is an 'object'? Each of these specifies some part of an object
		
		vector<Reactor *>reactorList; /* 
								* Determines part of ObjectSpec, 
							    * returns the products of a reaction between objects of type ObjectSpec. 
							    * Each Reactor encodes a particular type of reaction, and can be chained together in various ways
							    * Colliders feed the Reactors with compounds to try to react together (synchronously? asynchronously?), and the reactors
							    * decide what happens (and push the results to the world)
							    */
							   		
		vector<Collider *>colliderList; /*
								 * Depends on Arranger, returns a list of potentially reacting objects. 
								 * Each Collider can pipe to specific Reactors		
								 */
								 
		Arranger *spatialArranger; /*
						 * Spatial model - continuous, singular grid, 'bucket' grid, arbitrary network, etc. Determines ObjectSpec
						 * 
						 * The Arranger determines the relationship between spatialContainers, which are generic storage objects for compounds
						 * The spatialContainers can also hold certain optimized forms of their contents (e.g. hashes, etc), dependent on the
						 * details of the Reactor and Arranger.
						 */
		
		vector<Dynamic *>dynamicList; /*
								 * Different dynamical events that occur to redistribute things in the spatial model
								 * Each Dynamics module performs some re-arrangement of objects between spatialContainers
								 * Dynamics are dependent on the choice of Arranger
								 */
		
		vector<Initializer*>initializerList; // Sets up initial conditions. This can run per-site or per-particle or whatever.		
		vector<Analyzer*> analyzerList; // things that analyze the results
		vector<void*> Drivers; // things that mess with the simulation
		
		vector<void*> spatialContainers; // vector of distinct containers for objects
		vector<void*> objects; // Global list of objects
		
		vector<Interaction> getCollisionList();
		Simulation();
		Simulation(string filename);
};

vector<Interaction> Simulation::getCollisionList() // Get the collisions from the collider, given current data
{
	return Collider->getCollisionList();
}

Simulation::Simulation()
{
}

Simulation::Simulation(string filename)
{
}

int main(int argc, char **argv)
{
}
